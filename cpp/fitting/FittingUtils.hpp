//!
//! \file FittingUtils.hpp
//! \author Marco Casagrande
//! \date 28.06.2022
//!
//! \brief Utility class including algorithms used in- or after fitting.
//!

#pragma once

#include <cmath> 
#include <utility>

namespace fitting {

//!
//! \brief Utility class including algorithms used in- or after fitting.
//!
class FittingUtils
{
public:
    //!
    //! \brief Utility function to check if two points with floating point coordinates are approximately equal.
    //!
    //! \tparam TFloat  Type of floating point values.
    //!
    //! \param[in] pointA  Point A.
    //! \param[in] pointB  Point B.
    //!
    //! \return True, if point A and B are approximately equal.
    //!
    template<typename TFloat>
    static bool arePointsEqual(const std::pair<TFloat, TFloat>& pointA, const std::pair<TFloat, TFloat>& pointB);

    //!
    //! \brief Compute the point projection on the line defined by the slope and the offset.
    //!
    //! \tparam TFloat  Type of floating point values.
    //!
    //! \param[in] point  Point.
    //! \param[in] slope  Slope of the line.
    //! \param[in] offset  Offset of the line; if slope is infinte, offset should be in x direction.
    //!
    //! \return Projection of the point on the line.
    //!
    template<typename TFloat>
    static std::pair<TFloat, TFloat> computePointProjectionOnLine(const std::pair<TFloat, TFloat>& point,
                                                                  const std::pair<TFloat, TFloat>& pointA,
                                                                  const std::pair<TFloat, TFloat>& pointB);

    template<typename TFloat>
    static std::pair<TFloat, TFloat>
    computePointProjectionOnLine(const std::pair<TFloat, TFloat>& point, TFloat slope, TFloat offset);

    //!
    //! \brief Compute the slope of a line defined by two points.
    //!
    //! \tparam TFloat  Type of floating point values.
    //!
    //! \param[in] pointA  Point A.
    //! \param[in] pointB  Point B.
    //!
    //! \return Slope of the line; nan if the two points are the same and infinity if the two points share the same x coordinate.
    //!
    template<typename TFloat>
    static TFloat computeSlope(const std::pair<TFloat, TFloat>& pointA, const std::pair<TFloat, TFloat>& pointB);

    //!
    //! \brief Compute the squared distance between two points A and B.
    //!
    //! \tparam TFloat  Type of floating point values.
    //!
    //! \param[in] pointA  Point A.
    //! \param[in] pointB  Point B.
    //!
    //! \return Squared distance between two points.
    //!
    template<typename TFloat>
    static TFloat computeSquaredDistanceToPoint(const std::pair<TFloat, TFloat>& pointA,
                                                const std::pair<TFloat, TFloat>& pointB);

    //!
    //! \brief Compute the squared distance between a point and a segment defined by two points A and B.
    //!
    //! \tparam TFloat  Type of floating point values.
    //!
    //! \param[in] point  Point.
    //! \param[in] pointA  Point A.
    //! \param[in] pointB  Point B.
    //!
    //! \return Squared distance between the point and the segment.
    //!
    template<typename TFloat>
    static TFloat computeSquaredDistanceToSegment(const std::pair<TFloat, TFloat>& point,
                                                  const std::pair<TFloat, TFloat>& pointA,
                                                  const std::pair<TFloat, TFloat>& pointB);

    //!
    //! \brief Simplifyes a polyline defined by a range via the Ramer–Douglas–Peucker algorithm  by discarding the points that are closer than an epsilon value to a polyline segment.
    //! \details Based on the <a href="https://karthaus.nl/rdp/"> original implementation</a> and on its <a href="https://en.wikipedia.org/wiki/Ramer%E2%80%93Douglas%E2%80%93Peucker_algorithm#Pseudocode">pseudocode</a>.
    //!
    //! \tparam TFloat  Type of floating point values.
    //! \tparam Iterator  Point iterator.
    //! \tparam Get2DCoords  Functor that return a std::pair with the 2D coordinates of the points.
    //!
    //! \param[in] first  First point of the polyline.
    //! \param[in] last  Last point of the polyline.
    //! \param[in] get2DCoords  Functor that return a std::pair with the 2D coordinates of the points.
    //! \param[in] epsilon  Minimum distance between points and the closest polyline segment below which points are discarded; if less or equal to zero, no simplification is carried out.
    //!
    //! \return Iterator to the first discarded point; the relative order of the elements on the left side of the iterator is preserved while there is no guarantee about the order of the discarded points.
    //!
    template<typename TFloat, typename Iterator, typename Get2DCoords>
    static Iterator simplifyPolyline(Iterator first, Iterator last, Get2DCoords get2DCoords, TFloat epsilon);

    //!
    //! \brief Check if two floating point numbers are approximately equal within a set tolerance.
    //!
    //! \tparam TFloat  Type of floating point values.
    //!
    //! \param[in] floatA  First floating point number.
    //! \param[in] floatB  Second floating point number.
    //! \param[in] eps  Tolerance.
    //!
    //! \return True, if the absolute difference between the two floating point numbers is less than the tolerance.
    //!
    template<typename TFloat>
    static inline TFloat fuzzyCompare(TFloat floatA, TFloat floatB, TFloat eps = 1e-9)
    {
        return std::abs(floatA - floatB) < eps;
    }

};

} // namespace fitting

#include "FittingUtilsImpl.hpp"
