//!
//! \file FittingUtilsImpl.hpp
//! \author Marco Casagrande
//! \date 28.06.2022
//!

#pragma once

#include <algorithm>
#include <limits>

namespace fitting {

template<typename TFloat>
bool FittingUtils::arePointsEqual(const std::pair<TFloat, TFloat>& pointA, const std::pair<TFloat, TFloat>& pointB)
{
    return fuzzyCompare(pointA.first, pointB.first) && fuzzyCompare(pointA.second, pointB.second);
}

template<typename TFloat>
std::pair<TFloat, TFloat> FittingUtils::computePointProjectionOnLine(const std::pair<TFloat, TFloat>& point,
                                                                     const std::pair<TFloat, TFloat>& pointA,
                                                                     const std::pair<TFloat, TFloat>& pointB)
{
    if (arePointsEqual(pointA, pointB))
    {
        return std::make_pair(std::numeric_limits<TFloat>::quiet_NaN(), std::numeric_limits<TFloat>::quiet_NaN());
    }
    if (arePointsEqual(point, pointA))
    {
        return pointA;
    }
    if (arePointsEqual(point, pointB))
    {
        return pointB;
    }

    const auto slope{computeSlope(pointA, pointB)};

    const auto offset{(std::isinf(slope)) ? TFloat{pointA.first}
                                          : TFloat{pointA.second} - TFloat{pointA.first} * slope};

    if (fuzzyCompare(slope, TFloat{0.}))
    {
        return std::make_pair(point.first, offset);
    }
    if (std::isinf(slope))
    {
        return std::make_pair(offset, point.second);
    }

    const auto counterSlope{TFloat{-1.} / slope};
    const auto counterOffset{TFloat{point.second} - TFloat{point.first} * counterSlope};
    const auto x{(offset - counterOffset) / (counterSlope - slope)};
    const auto y{counterSlope * x + counterOffset};
    return std::make_pair(x, y);
}

template<typename TFloat>
std::pair<TFloat, TFloat>
FittingUtils::computePointProjectionOnLine(const std::pair<TFloat, TFloat>& point, TFloat slope, TFloat offset)
{
    if (fuzzyCompare(slope, TFloat{0.}))
    {
        return std::make_pair(point.first, offset);
    }
    if (std::isinf(slope))
    {
        return std::make_pair(offset, point.second);
    }
    const auto counterSlope{TFloat{-1.} / slope};
    const auto counterOffset{point.second - point.first * counterSlope};
    const auto x{(offset - counterOffset) / (counterSlope - slope)};
    const auto y{counterSlope * x + counterOffset};
    return std::make_pair(x, y);
}

template<typename TFloat>
TFloat FittingUtils::computeSlope(const std::pair<TFloat, TFloat>& pointA, const std::pair<TFloat, TFloat>& pointB)
{
    return (arePointsEqual(pointA, pointB)) ? std::numeric_limits<TFloat>::quiet_NaN()
                                            : (pointB.second - pointA.second) / (pointB.first - pointA.first);
}

template<typename TFloat>
TFloat FittingUtils::computeSquaredDistanceToPoint(const std::pair<TFloat, TFloat>& pointA,
                                                   const std::pair<TFloat, TFloat>& pointB)
{
    return std::pow(pointA.first - pointB.first, TFloat{2.}) + std::pow(pointA.second - pointB.second, TFloat{2.});
}

template<typename TFloat>
TFloat FittingUtils::computeSquaredDistanceToSegment(const std::pair<TFloat, TFloat>& point,
                                                     const std::pair<TFloat, TFloat>& pointA,
                                                     const std::pair<TFloat, TFloat>& pointB)
{
    const auto slope{computeSlope(pointA, pointB)};

    const auto offset{(std::isinf(slope)) ? TFloat{pointA.first}
                                          : TFloat{pointA.second} - TFloat{pointA.first} * slope};

    const auto pointProjected{computePointProjectionOnLine(point, slope, offset)};
    const auto squaredDistanceToLine{computeSquaredDistanceToPoint(point, pointProjected)};
    return (pointProjected.first > std::min(pointA.first, pointB.first)
            && pointProjected.first < std::max(pointA.first, pointB.first))
               ? squaredDistanceToLine
               : std::min(computeSquaredDistanceToPoint(point, pointA), computeSquaredDistanceToPoint(point, pointB));
}

template<typename TFloat, typename Iterator, typename Get2DCoords>
Iterator FittingUtils::simplifyPolyline(Iterator first, Iterator last, Get2DCoords get2DCoords, TFloat epsilon)
{
    auto pointA{first};
    auto pointB{std::prev(last)};

    // end recursion if two or less points are available, if start and and points are the same or epsilion is not positive
    if (std::distance(pointA, pointB) < 2 || *pointA == *pointB || epsilon <= TFloat{0})
    {
        return last;
    }

    const auto pointFurthest{std::max_element(std::next(pointA), pointB, [&](const auto& pt1, const auto& pt2) {
        return computeSquaredDistanceToSegment<TFloat>(get2DCoords(pt1), get2DCoords(*pointA), get2DCoords(*pointB))
               < computeSquaredDistanceToSegment<TFloat>(get2DCoords(pt2), get2DCoords(*pointA), get2DCoords(*pointB));
    })};

    if (computeSquaredDistanceToSegment<TFloat>(get2DCoords(*pointFurthest), get2DCoords(*pointA), get2DCoords(*pointB))
        > std::pow(epsilon, TFloat{2}))
    {
        // if max distance is greater than epsilon, recursively simplify left and right half of the polyline
        // important to start from right not to swap iterator to max before needed
        pointB = simplifyPolyline<TFloat>(pointFurthest, std::next(pointB), get2DCoords, epsilon);
        pointA = simplifyPolyline<TFloat>(pointA, std::next(pointFurthest), get2DCoords, epsilon);

        // partition points so that the ones to keep are moved on the left side of the retutn iterator preserving their relative order
        // while the ones to discard are on the right
        auto pt{std::next(pointFurthest)};
        while (pt != pointB)
        {
            std::iter_swap(pointA++, pt++);
        }
        return pointA;
    }
    else
    {
        // else, keep only first and last point by swapping the latter to the second position
        std::iter_swap(pointB, std::next(pointA));
        return std::next(pointA, 2);
    }
}

} // namespace fitting
