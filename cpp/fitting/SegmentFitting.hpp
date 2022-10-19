//!
//! \file SegmentFitting.hpp
//! \author Marco Casagrande
//! \date 05.05.2022
//!
//! \brief Fit 2D segments to sets of points.
//!

#pragma once

#include <utility>

namespace fitting {

//!
//! \brief Segment 2D.
//!
//! \tparam TFloat  Type of floating point values.
//!
template<typename TFloat>
struct Segment2D
{
    //!
    //! \brief Default constructor.
    //!
    Segment2D() = default;

    //!
    //! \brief Construct a new Segment2D object.
    //!
    //! \param[in] begin  2D coordinates of the start of the segment.
    //! \param[in] end  2D coordinates of the end of the segment.
    //! \param[in] mse  Mean-Squared-Error of the estimate.
    //!
    Segment2D(const std::pair<TFloat, TFloat>& begin, const std::pair<TFloat, TFloat>& end, TFloat mse)
      : m_begin{begin}, m_end{end}, m_mse{mse}, m_isValid{true} {};

    //! 2D coordinates of the start of the segment.
    std::pair<TFloat, TFloat> m_begin;

    //! 2D coordinates of the end of the segment.
    std::pair<TFloat, TFloat> m_end;

    //! Mean-Squared-Error of the estimate.
    TFloat m_mse;

    //! Validity flag.
    bool m_isValid;
};

//!
//! \brief class for the clustering of point by a DBSCAN algorithm
//!
class SegmentFitting
{
public:
    //!
    //! \brief Fit segment 2D to points using the RANSAC approach.
    //!
    //! \tparam TFloat  Type of floating point values.
    //! \tparam Iterator  Type of iterator.
    //! \tparam Get2DCoords  Functor that return a std::pair with the 2D coordinates of the points.
    //!
    //! \param[in] first  Iterator to first point in the range to fit the segment to.
    //! \param[in] last  Iterator to last point in the range to fit the segment to.
    //! \param[in] get2DCoords  Getter functor instance.
    //! \param[in] nIterations  Number of RANSAC iterations; should be > 0.
    //! \param[in] nSamples  Number of points to sample to validate the estimate; use all points by default.
    //! \param[in] maxInliersDistance  Maximum inlier distance from segment above which outliers are discarded; by default all sampled points are considered inliers.
    //!
    //! \return Segment2D fit to the points.
    //!
    template<typename TFloat, typename Iterator, typename Get2DCoords>
    static Segment2D<TFloat> fitSegment2D(Iterator first,
                                          Iterator last,
                                          Get2DCoords get2DCoords,
                                          size_t nIterations,
                                          size_t nSamples           = 0,
                                          TFloat maxInliersDistance = TFloat{0.});

    //!
    //! \brief Partition 2D points at a set maximum distance from the segment, by moving them on the left side of the range.
    //!
    //! \tparam TFloat  Type of floating point values.
    //! \tparam Iterator  Type of non-const iterator.
    //! \tparam Get2DCoords  Functor that return a std::pair with the 2D coordinates of the points.
    //!
    //! \param[in] first  Non-const iterator to first point in the range to partition.
    //! \param[in] last  Non-const iterator to last point in the range to partition.
    //! \param[in] get2DCoords  Getter functor instance.
    //! \param[in] segment  Segment to check for inliers.
    //! \param[in] maxInliersDistance  Maximum distance from segment of inliers; if set to 0, all points are considered inliers.
    //!
    //! \return Iterator to the first outlier.
    //!
    template<typename TFloat, typename Iterator, typename Get2DCoords>
    static Iterator partitionInliers2D(Iterator first,
                                       Iterator last,
                                       Get2DCoords get2DCoords,
                                       const Segment2D<TFloat>& segment,
                                       TFloat maxInliersDistance);
};

} // namespace fitting

#include "SegmentFittingImpl.hpp"
