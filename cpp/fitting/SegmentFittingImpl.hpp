//!
//! \file SegmentFittingImpl.hpp
//! \author Marco Casagrande
//! \date 05.05.2022
//!

#pragma once

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <random>

namespace fitting {

template<typename TFloat, typename Iterator, typename Get2DCoords>
Segment2D<TFloat> SegmentFitting::fitSegment2D(Iterator first,
                                               Iterator last,
                                               Get2DCoords get2DCoords,
                                               size_t nIterations,
                                               size_t nSamples,
                                               TFloat maxInliersDistance)
{
    const auto nPoints{static_cast<size_t>(std::distance(first, last))};
    if (nPoints < 2)
    {
        logWarning("Points are not enough to fit segment!");
        return Segment2D<TFloat>{};
    }
    if (nIterations == 0)
    {
        logWarning("Number of iterations is 0! Defaulting to 1... Estimate is likely inaccurate!");
        nIterations = 1;
    }

    nSamples = (nSamples == 0) ? nPoints : std::min(nSamples, nPoints);

    std::random_device rd;
    std::mt19937 gen{rd()};

    const auto maxSquaredDistance{(maxInliersDistance != TFloat{0.}) ? std::pow(maxInliersDistance, TFloat{2.})
                                                                     : std::numeric_limits<TFloat>::infinity()};
    Segment2D<TFloat> bestFitSegment;
    auto smallestSumSquaredDistance{std::numeric_limits<TFloat>::infinity()};

    while (nIterations-- > 0)
    {
        const auto iA{std::uniform_int_distribution<>{0ul, nPoints - 1ul}(gen)};
        size_t iB;
        do
        {
            iB = std::uniform_int_distribution<>{0ul, nPoints - 1ul}(gen);
        } while (iA == iB);

        const auto pointA{get2DCoords(*(first + iA))};
        const auto pointB{get2DCoords(*(first + iB))};

        const auto slope{FittingUtils::computeSlope(pointA, pointB)};
        if (std::isnan(FittingUtils::computeSlope(pointA, pointB)))
        {
            continue;
        }

        // use x as offset for distance computation fr longitudinal lines
        const auto offset{(std::isinf(slope)) ? pointA.first : pointA.second - pointA.first * slope};

        std::vector<std::pair<TFloat, TFloat>> segmentPoints{pointA, pointB};
        segmentPoints.reserve(nSamples + 2ul);

        TFloat sumSquaredDistance{0.};

        // samples without replacement using Robert Floyd algorithm
        auto n{nPoints - nSamples};
        std::vector<bool> wasSampled(nPoints);

        while (n < nPoints && sumSquaredDistance < smallestSumSquaredDistance)
        {
            auto i{std::uniform_int_distribution<>{0ul, n}(gen)};
            i             = (wasSampled.at(i)) ? n : i;
            wasSampled[i] = true;
            n++;

            const auto point{get2DCoords(*(first + i))};
            auto pointOnSegment{FittingUtils::computePointProjectionOnLine(point, slope, offset)};
            const auto squaredDistance{FittingUtils::computeSquaredDistanceToPoint(point, pointOnSegment)};

            sumSquaredDistance += std::min(squaredDistance, maxSquaredDistance);

            if (squaredDistance <= maxSquaredDistance)
            {
                segmentPoints.emplace_back(std::move(pointOnSegment));
            }
        }

        if (sumSquaredDistance >= smallestSumSquaredDistance)
        {
            continue;
        }

        smallestSumSquaredDistance = sumSquaredDistance;

        const auto segmentBeginEnd{
            std::minmax_element(segmentPoints.begin(), segmentPoints.end(), [](const auto& p1, const auto& p2) {
                return p1.first < p2.first;
            })};

        bestFitSegment
            = Segment2D<TFloat>{*segmentBeginEnd.first, *segmentBeginEnd.second, sumSquaredDistance / nSamples};
    }

    return bestFitSegment;
}

template<typename TFloat, typename Iterator, typename Get2DCoords>
Iterator SegmentFitting::partitionInliers2D(Iterator first,
                                            Iterator last,
                                            Get2DCoords get2DCoords,
                                            const Segment2D<TFloat>& segment,
                                            TFloat maxInliersDistance)
{
    if (maxInliersDistance == TFloat{0.})
    {
        return last;
    }

    if (std::isnan(FittingUtils::computeSlope(segment.m_begin, segment.m_end)))
    {
        return first;
    }

    const auto maxSquaredDistance{std::pow(maxInliersDistance, TFloat{2.})};
    return std::partition(first, last, [&](auto point) {
        return FittingUtils::computeSquaredDistanceToSegment(get2DCoords(point), segment.m_begin, segment.m_end)
               < maxSquaredDistance;
    });
}

} // namespace fitting
