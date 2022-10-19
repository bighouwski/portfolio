//!
//! \file Example.hpp
//! \author Marco Casagrande
//! \date 28.06.2022
//!
//! \brief Class to fit polylines to an image by tracing the topological skeleton on a vector representing 2D image data.
//! \details This C++ implementation of the skeleton tracing algorithm is based on the <a href="https://github.com/LingDong-/skeleton-tracing">original C implementation</a> from Lingdong Huang.
//!

#pragma once


#include <algorithm>
#include <iterator>
#include <vector>

namespace fitting {

//!
//! \brief Class to fit polylines to an image by tracing the topological skeleton on a vector representing 2D image data.
//! \details This C++ implementation of the skeleton tracing algorithm is based on the <a href="https://github.com/LingDong-/skeleton-tracing">original C implementation</a> from Lingdong Huang.
//!
class PolylineFitting
{
public:
    //!
    //! \brief Fits polylines to 2D image data based on the provided unary predicate provided to determine the "on"/"off" state of the pixels.
    //!
    //! \tparam TP  Type of pixels.
    //! \tparam UnaryPredicate  Returns true if pixel is to be considered "on".
    //!
    //! \param[in] imageData  Image data of custom type; if the "doThinning" option is selected, the image is also thinned before fitting the polylines.
    //! \param[in] rows  Number of rows of the provided image; the size of image data should be equal to the number of rows multipled by the number of columns.
    //! \param[in] cols  Number of coulumns of the provided image; the size of image data should be equal to the number of rows multipled by the number of columns.
    //! \param[in] isPixelOn  Unary predicate returning true if any given pixel is to be considered "on".
    //! \param[in] minSectionSize  Size of the smallest image chunks to fit polylines to; smaller chunks lead to a higher resolution and to potentially more noisy polylines.
    //! \param[in] maxRecursions  Max number of times the algorithm should keep recursively splitting the image; if zero, keep recursing until minimum section size is reached or no more "on" pixel are found.
    //! \param[in] doThinning  Thins the original image before fitting the polylines; not necessary if image has thin pixel strokes already.
    //!
    //! \return Polylines fitted to the image.
    //!
    template<typename TP, typename UnaryPredicate>
    static std::vector<std::vector<std::pair<int32_t, int32_t>>> fitPolylines(const TP* imageData,
                                                                              size_t rows,
                                                                              size_t cols,
                                                                              UnaryPredicate isPixelOn,
                                                                              size_t minSectionSize = 3,
                                                                              size_t maxRecursions  = 0,
                                                                              bool doThinning       = true)
    {
        return fitPolylines(BitImage{imageData, isPixelOn, static_cast<int32_t>(rows), static_cast<int32_t>(cols)},
                            minSectionSize,
                            maxRecursions,
                            doThinning);
    }

private:
    //!
    //! \brief Helper class to handle 2D binary image data.
    //!
    class BitImage : public std::vector<uint8_t>
    {
    public:
        //!
        //! \brief Construct a new BitImage object with the specified size.
        //!
        //! \param[in] rows  Number of rows.
        //! \param[in] cols  Number of coulumns.
        //!
        BitImage(int32_t rows, int32_t cols)
          : std::vector<uint8_t>(static_cast<size_t>(rows * cols)), m_rows{rows}, m_cols{cols} {};

        //!
        //! \brief Construct a new BitImage object from existing image data evaluated according to the provided unary predicate.
        //!
        //! \tparam TP  Type of pixels.
        //! \tparam UnaryPredicate  Returns true if pixel is to be considered "on".
        //!
        //! \param[in] imageData  Image data of custom type.
        //! \param[in] isPixelOn  Unary predicate returning true if any given pixel is to be considered "on".
        //! \param[in] rows  Number of rows of the provided image; the size of image data needs to be equal to the number of rows multipled by the number of columns.
        //! \param[in] cols  Number of coulumns of the provided image; the size of image data needs to be equal to the number of rows multipled by the number of columns.
        //!
        template<typename TP, typename UnaryPredicate>
        BitImage(const TP* imageData, UnaryPredicate isPixelOn, int32_t rows, int32_t cols) : BitImage{rows, cols}
        {
            std::transform(imageData, imageData + rows * cols, begin(), isPixelOn);
        };

        //!
        //! \brief Retrieve the 2D coordinates of the provided pixel iterator.
        //!
        //! \param[in] px  Iterator to pixel.
        //!
        //! \return 2D coordinates of the pixel.
        //!
        inline std::pair<int32_t, int32_t> coords(const_iterator px) const noexcept
        {
            return std::make_pair(std::distance(begin(), px) / m_cols, std::distance(begin(), px) % m_cols);
        };

        //!
        //! \brief Retrieves the const iterator to the target pixel located at the coordinates provided.
        //!
        //! \param[in] row  Row of the target pixel.
        //! \param[in] col  Column of the target pixel.
        //!
        //! \return Const iterator to the target pixel.
        //!
        inline const_iterator find(int32_t row, int32_t col) const noexcept { return begin() + m_cols * row + col; };

        //!
        //! \brief Retrieves the iterator to the target pixel located at the coordinates provided.
        //!
        //! \param[in] row  Row of the target pixel.
        //! \param[in] col  Column of the target pixel.
        //!
        //! \return Const iterator to the target pixel.
        //!
        inline iterator find(int32_t row, int32_t col) noexcept { return begin() + m_cols * row + col; };

        //!
        //! \brief Retrieves the const iterators of the pixel within a subsection of the image.
        //!
        //! \param[in] r0  Starting row of the image setion.
        //! \param[in] c0  Starting column of the image section.
        //! \param[in] rows  Number of rows of the image section; if negative, retrieve backwards.
        //! \param[in] cols  Number of columns of the image section; if negative, retrieve backwards.
        //!
        //! \return Const iterators to the pixels within the image section.
        //!
        inline std::vector<const_iterator> section(int32_t r0, int32_t c0, int32_t rows, int32_t cols) const noexcept
        {
            std::vector<const_iterator> pixels;
            pixels.reserve(static_cast<size_t>(rows * cols));
            for (auto r{r0}; r < r0 + rows; ++r)
            {
                for (auto c{c0}; c < c0 + cols; ++c)
                {
                    pixels.emplace_back(begin() + m_cols * r + c);
                }
            }
            return pixels;
        };

        //!
        //! \brief Retrieves the iterators of the pixel within a subsection of the image.
        //!
        //! \param[in] r0  Starting row of the image setion.
        //! \param[in] c0  Starting column of the image section.
        //! \param[in] rows  Number of rows of the image section; if negative, retrieve backwards.
        //! \param[in] cols  Number of columns of the image section; if negative, retrieve backwards.
        //!
        //! \return Iterators to the pixels within the image section.
        //!
        inline std::vector<iterator> section(int32_t r0, int32_t c0, int32_t rows, int32_t cols) noexcept
        {
            std::vector<iterator> pixels;
            pixels.reserve(static_cast<size_t>(rows * cols));
            for (auto r{r0}; r < r0 + rows; ++r)
            {
                for (auto c{c0}; c < c0 + cols; ++c)
                {
                    pixels.emplace_back(begin() + m_cols * r + c);
                }
            }
            return pixels;
        };

    public:
        //! Number of image rows.
        const int32_t m_rows;

        //! Number of image columns.
        const int32_t m_cols;
    };

private:
    //!
    //! \brief Entry point that initiates the recursive algorithm.
    //!
    //! \param[in] image  BitImage to trace the skeleton of.
    //! \param[in] minSectionSize  Size of the smallest image chunks to fit polylines to; smaller chunks lead to a higher resolution and to potentially more noisy polylines.
    //! \param[in] maxRecursions  Max number of times the algorithm should keep recursively splitting the image; if zero, keep recursing until minimum section size is reached or no more "on" pixel are found.
    //! \param[in] doThinning  Thins the original image before fitting the polylines; not necessary if image has thin pixel strokes already.
    //!
    //! \return Polylines fit to the image section.
    //!
    static std::vector<std::vector<std::pair<int32_t, int32_t>>>
    fitPolylines(BitImage&& image, size_t minSectionSize, size_t maxRecursions, bool doThinning);

    //!
    //! \brief Main algorithm to fit the polylines to a 2D image section by tracing its topological skeleton.
    //! \details The algorithm performs the following steps recursively:
    //! - splits the image into increasingly smaller sections on the row/column with the least "on" pixels and closest to the center;
    //! - fits segments to the sections as their size becomes small enough or the number of iteration gets larger than a maximum threshold;
    //! - merges segments into increasingly longer polylines until the top of the call stack is reached.
    //!
    //! \param[in] image  BitImage to trace the skeleton of.
    //! \param[in] r0  Starting row of the image setion.
    //! \param[in] c0  Starting column of the image section.
    //! \param[in] rows  Number of rows of the image section.
    //! \param[in] cols  Number of columns of the image section.
    //! \param[in] iter  Current iteration.
    //! \param[in] minSectionSize  Minimum size of the image sections to which the polyline segments are fit to.
    //! \param[in] maxRecursions  Maximum number of iterations in which the image can be further split into sections.
    //!
    //! \return Polylines fit to the image section.
    //!
    static std::vector<std::vector<BitImage::const_iterator>> fitPolylines(const BitImage& image,
                                                                           int32_t r0,
                                                                           int32_t c0,
                                                                           int32_t rows,
                                                                           int32_t cols,
                                                                           int32_t iter,
                                                                           int32_t minSectionSize,
                                                                           int32_t maxRecursions);

    //!
    //! \brief Fit segments to a 2D image section from its frame pixels to the most likely intersection pixel.
    //!
    //! \param[in] image  BitImage to fit the segments to.
    //! \param[in] r0  Starting row of the image setion.
    //! \param[in] c0  Starting column of the image section.
    //! \param[in] rows  Number of rows of the image section.
    //! \param[in] cols  Number of columns of the image section.
    //!
    //! \return Segments fitted to the image section
    //!
    static std::vector<std::vector<BitImage::const_iterator>>
    fitSegments(const BitImage& image, int32_t r0, int32_t c0, int32_t rows, int32_t cols);

    //!
    //! \brief Merge polylines on their common junction pixel.
    //!
    //! \param[in,out] destPolylines  Destination polylines to merge the source polylines to.
    //! \param[in] srcPolylines  Source polylines to merge into the destination polylines.
    //!
    static std::vector<std::vector<BitImage::const_iterator>>
    mergePolylines(std::vector<std::vector<BitImage::const_iterator>>&& destPolylines,
                   std::vector<std::vector<BitImage::const_iterator>>&& srcPolylines);

    //!
    //! \brief Thins the 2D image in-place, so that the image is reduced to its toplogical skeleton with a 1-pixel thickness.
    //! \details This C++ implementation is based on the <a href="http://agcggs680.pbworks.com/f/Zhan-Suen_algorithm.pdf>original algorithm</a> from T. Y. Zhang and C. Y. Suen.
    //!
    //! \param[in] image  BitImage to be thinned.
    //!
    static void thinImage(BitImage& image);
};

} // namespace fitting
