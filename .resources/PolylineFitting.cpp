//!
//! \file PolylineFitting.cpp
//! \author Marco Casagrande
//! \date 28.06.2022
//!

#include "PolylineFitting.hpp"

#include <cmath>
#include <numeric>

namespace fitting {

std::vector<std::vector<std::pair<int32_t, int32_t>>>
PolylineFitting::fitPolylines(BitImage&& image, size_t minSectionSize, size_t maxRecursions, bool doThinning)
{
    const auto minSize{3ul};
    if (static_cast<size_t>(image.m_rows) < minSize || static_cast<size_t>(image.m_cols) < minSize)
    {
        logWarning("Impossible to fit polylines to an image smaller than minimum size of 3x3!");
        return std::vector<std::vector<std::pair<int32_t, int32_t>>>{};
    }

    if (doThinning)
    {
        thinImage(image);
    }

    const auto pxPolylines{
        fitPolylines(image,
                     0,
                     0,
                     image.m_rows,
                     image.m_cols,
                     0,
                     static_cast<int32_t>(std::max(minSectionSize, minSize)),
                     maxRecursions != 0ul ? static_cast<int32_t>(maxRecursions) : std::numeric_limits<int32_t>::max())};

    std::vector<std::vector<std::pair<int32_t, int32_t>>> coordsPolylines;
    for (const auto& pl : pxPolylines)
    {
        coordsPolylines.emplace_back(pl.size());
        std::transform(pl.begin(), pl.end(), coordsPolylines.back().begin(), [&](auto px) { return image.coords(px); });
    }
    return coordsPolylines;
}

std::vector<std::vector<PolylineFitting::BitImage::const_iterator>>
PolylineFitting::fitPolylines(const PolylineFitting::BitImage& image,
                              int32_t r0,
                              int32_t c0,
                              int32_t rows,
                              int32_t cols,
                              int32_t iter,
                              int32_t minSectionSize,
                              int32_t maxRecursions)
{
    // lambda to check if image section has at least one "on" pixel
    auto isEmpty{
        [](const auto& pixels) { return !std::any_of(pixels.begin(), pixels.end(), [](auto px) { return *px; }); }};

    // end recursion if image section does not contain any "on" pixels
    if (isEmpty(image.section(r0, c0, rows, cols)))
    {
        return std::vector<std::vector<PolylineFitting::BitImage::const_iterator>>{};
    }

    // split only makes sense with at least 5 rows or columns, so that the polylines can be merged at the center
    const auto minSplitSize{std::max(minSectionSize, 5)};

    // end recursion when image section is smaller than minimum size or number of iterations exceeds the maximum threshold
    if (iter >= maxRecursions || (cols < minSplitSize && rows < minSplitSize))
    {
        return fitSegments(image, r0, c0, rows, cols);
    }

    int32_t rSplit{-1};
    int32_t cSplit{-1};
    auto& bestSplit{rows >= cols ? rSplit : cSplit};
    auto minOn{std::numeric_limits<int64_t>::max()};
    const auto iMax{rows >= cols ? rows - 4 : cols - 4};

    // identify best splitting candidate which is closer to center and with the least "on" pixel
    for (auto i{0}; i < iMax; i++)
    {
        // increasingly larger offset with alternating sign at each iteration, e.g. 0, -1, 1, -2, 2, ...
        const auto offset{static_cast<int32_t>(std::pow(-1, i)) * (i + 1) / 2};
        const auto splitCandidate{rows >= cols ? r0 + rows / 2 + offset : c0 + cols / 2 + offset};

        const auto pixels{rows >= cols ? image.section(splitCandidate, c0, 1, cols)
                                       : image.section(r0, splitCandidate, rows, 1)};

        const auto nOnPixels{std::count_if(pixels.begin(), pixels.end(), [](const auto& px) { return *px; })};

        bestSplit = nOnPixels < minOn ? splitCandidate : bestSplit;
        minOn     = std::min(minOn, nOnPixels);

        // early stopping if no "on" pixel
        if (minOn == 0)
        {
            break;
        }
    }

    // recursively split image section in two and fit polylines to each subsection
    if (bestSplit == rSplit)
    {
        return mergePolylines(
            fitPolylines(image, r0, c0, rSplit - r0 + 1, cols, iter + 1, minSectionSize, maxRecursions),
            fitPolylines(image, rSplit, c0, r0 + rows - rSplit, cols, iter + 1, minSectionSize, maxRecursions));
    }
    else // bestSplit == cSplit
    {
        return mergePolylines(
            fitPolylines(image, r0, c0, rows, cSplit - c0 + 1, iter + 1, minSectionSize, maxRecursions),
            fitPolylines(image, r0, cSplit, rows, c0 + cols - cSplit, iter + 1, minSectionSize, maxRecursions));
    }
}

std::vector<std::vector<PolylineFitting::BitImage::const_iterator>>
PolylineFitting::fitSegments(const PolylineFitting::BitImage& image, int32_t r0, int32_t c0, int32_t rows, int32_t cols)
{
    // bottom-right image section row and column
    const auto r1{r0 + rows - 1};
    const auto c1{c0 + cols - 1};

    // retrieve image section frame clockwise
    auto topRow{image.section(r0, c0, 1, cols - 1)};
    auto rightCol{image.section(r0, c1, rows - 1, 1)};
    auto bottomRow{image.section(r1, c0 + 1, 1, cols - 1)};
    auto leftCol{image.section(r0 + 1, c0, rows - 1, 1)};

    std::vector<PolylineFitting::BitImage::const_iterator> pixelsFrame;
    pixelsFrame.reserve(static_cast<size_t>(2 * (rows + cols)));
    pixelsFrame.insert(
        pixelsFrame.end(), std::make_move_iterator(topRow.begin()), std::make_move_iterator(topRow.end()));
    pixelsFrame.insert(
        pixelsFrame.end(), std::make_move_iterator(rightCol.begin()), std::make_move_iterator(rightCol.end()));
    pixelsFrame.insert(
        pixelsFrame.end(), std::make_move_iterator(bottomRow.rbegin()), std::make_move_iterator(bottomRow.rend()));
    pixelsFrame.insert(
        pixelsFrame.end(), std::make_move_iterator(leftCol.rbegin()), std::make_move_iterator(leftCol.rend()));

    auto pxOff{std::find_if(pixelsFrame.begin(), pixelsFrame.end(), [](auto px) { return *px ? false : true; })};
    auto pxOn{std::find_if(pixelsFrame.begin(), pixelsFrame.end(), [](auto px) { return *px ? true : false; })};

    // impossible to determine segments if all pixels are either "on" or "off"
    if (pxOff == pixelsFrame.end() || pxOn == pixelsFrame.end())
    {
        return std::vector<std::vector<PolylineFitting::BitImage::const_iterator>>{};
    }

    // rotate frame so that it starts from an "off" pixel
    std::rotate(pixelsFrame.begin(), pxOff, pixelsFrame.end());

    const auto pxCenter{image.find(r0 + rows / 2, c0 + cols / 2)};

    // iterate through image frame and fit segments from middle of "on" pixels to center of image
    std::vector<std::vector<PolylineFitting::BitImage::const_iterator>> segments;
    do
    {
        pxOn  = std::find_if(pxOff, pixelsFrame.end(), [&](auto px) { return *px ? true : false; });
        pxOff = std::find_if(pxOn, pixelsFrame.end(), [&](auto px) { return *px ? false : true; });
        const auto pxMid{*(pxOn + static_cast<int32_t>(pxOff - pxOn) / 2)};

        if (pxOn != pixelsFrame.end())
        {
            segments.emplace_back(std::vector<PolylineFitting::BitImage::const_iterator>{pxMid, pxCenter});
        }
    } while (pxOn != pixelsFrame.end());

    // if two segments, merge them into a longer segment
    if (segments.size() == 2)
    {
        return std::vector<std::vector<PolylineFitting::BitImage::const_iterator>>{
            std::vector<PolylineFitting::BitImage::const_iterator>{segments.front().front(), segments.back().front()}};
    }

    // if more than 2 segments, guess the intersection
    auto pxsCanvas{image.section(r0 + 1, c0 + 1, rows - 2, cols - 2)};

    // sort pixels from closest to farthest from center
    std::sort(pxsCanvas.begin(), pxsCanvas.end(), [&](const auto& pxA, const auto& pxB) {
        const auto coordsA{image.coords(pxA)};
        const auto coordsB{image.coords(pxB)};
        const auto coordsC{image.coords(pxCenter)};

        return std::abs(coordsA.first - coordsC.first) + std::abs(coordsA.second - coordsC.second)
               < std::abs(coordsB.first - coordsC.first) + std::abs(coordsB.second - coordsC.second);
    });

    // heuristic to stop early if enough pixels are on
    const auto minOnPixels{5};
    auto maxConv{-1l};
    BitImage::const_iterator pxIntersection;
    for (auto px : pxsCanvas)
    {
        // use convolution to find 3x3 pxel neighbourhood with the most "on" pixels closest to the center
        const auto coords{image.coords(px)};
        const auto pixels{image.section(coords.first - 1, coords.second - 1, 3, 3)};
        const auto conv{std::count_if(pixels.begin(), pixels.end(), [](auto pxx) { return *pxx; })};
        pxIntersection = conv > maxConv ? px : pxIntersection;
        maxConv        = std::max(maxConv, conv);

        // early stopping if pixel is good enough candidate
        if (maxConv >= minOnPixels)
        {
            break;
        }
    }

    // update segments with new estimated intersection
    for (auto& polyline : segments)
    {
        polyline.back() = pxIntersection;
    }

    return segments;
}

std::vector<std::vector<PolylineFitting::BitImage::const_iterator>>
PolylineFitting::mergePolylines(std::vector<std::vector<PolylineFitting::BitImage::const_iterator>>&& destPolylines,
                                std::vector<std::vector<PolylineFitting::BitImage::const_iterator>>&& srcPolylines)
{
    if (destPolylines.empty())
    {
        return srcPolylines;
    }
    if (srcPolylines.empty())
    {
        return destPolylines;
    }

    auto srcEnd{srcPolylines.end()};

    // compare the extremities of each destination - and source polylines not yet merged and merge them if they are the same
    for (auto& dest : destPolylines)
    {
        auto srcIt{srcEnd};

        // check for front-front extremities
        srcIt = std::find_if(
            srcPolylines.begin(), srcEnd, [&dest](const auto& src) { return dest.front() == src.front(); });
        if (srcIt != srcEnd)
        {
            dest.insert(
                dest.begin(), std::make_move_iterator(srcIt->rbegin()), std::make_move_iterator(srcIt->rend() - 1));
            std::iter_swap(srcIt, --srcEnd);
            continue;
        }

        // check for front-back extremities
        srcIt = std::find_if(
            srcPolylines.begin(), srcEnd, [&dest](const auto& src) { return dest.front() == src.back(); });
        if (srcIt != srcEnd)
        {
            dest.insert(
                dest.begin(), std::make_move_iterator(srcIt->begin()), std::make_move_iterator(srcIt->end() - 1));
            std::iter_swap(srcIt, --srcEnd);
            continue;
        }

        // check for back-front extremities
        srcIt = std::find_if(
            srcPolylines.begin(), srcEnd, [&dest](const auto& src) { return dest.back() == src.front(); });
        if (srcIt != srcEnd)
        {
            dest.insert(dest.end(), std::make_move_iterator(srcIt->begin() + 1), std::make_move_iterator(srcIt->end()));
            std::iter_swap(srcIt, --srcEnd);
            continue;
        }

        // check for back-back extremities
        srcIt = std::find_if(
            srcPolylines.begin(), srcEnd, [&dest](const auto& src) { return dest.back() == src.back(); });
        if (srcIt != srcEnd)
        {
            dest.insert(
                dest.end(), std::make_move_iterator(srcIt->rbegin() + 1), std::make_move_iterator(srcIt->rend()));
            std::iter_swap(srcIt, --srcEnd);
            continue;
        }
    }

    // insert unmerged source polylines
    destPolylines.insert(
        destPolylines.end(), std::make_move_iterator(srcPolylines.begin()), std::make_move_iterator(srcEnd));

    return destPolylines;
}

void PolylineFitting::thinImage(PolylineFitting::BitImage& image)
{
    // iterate through the image canvas and delete the pixels according to the original algorithm cited in the header until nothing left to delete
    auto difference{true};
    auto iteration{false};

    auto onPixels{image.section(1, 1, image.m_rows - 2, image.m_cols - 2)};

    while (difference)
    {
        onPixels.erase(std::remove_if(onPixels.begin(), onPixels.end(), [](auto px) { return *px == 0; }),
                       onPixels.end());

        std::vector<PolylineFitting::BitImage::iterator> flaggedPixels;
        flaggedPixels.reserve(onPixels.size());

        for (auto px : onPixels)
        {
            int32_t r;
            int32_t c;
            std::tie(r, c) = image.coords(px);

            const auto p2{*image.find(r - 1, c + 0)};
            const auto p3{*image.find(r - 1, c + 1)};
            const auto p4{*image.find(r + 0, c + 1)};
            const auto p5{*image.find(r + 1, c + 1)};
            const auto p6{*image.find(r + 1, c + 0)};
            const auto p7{*image.find(r + 1, c - 1)};
            const auto p8{*image.find(r + 0, c - 1)};
            const auto p9{*image.find(r - 1, c - 1)};

            const auto A{int32_t{!p2 && p3} + int32_t{!p3 && p4} + int32_t{!p4 && p5} + int32_t{!p5 && p6}
                         + int32_t{!p6 && p7} + int32_t{!p7 && p8} + int32_t{!p8 && p9} + int32_t{!p9 && p2}};

            const auto B{p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9};
            const auto m1{iteration ? p2 && p4 && p8 : p2 && p4 && p6};
            const auto m2{iteration ? p2 && p6 && p8 : p4 && p6 && p8};

            if (A == 1 && (B >= 2 && B <= 6) && !m1 && !m2)
            {
                flaggedPixels.emplace_back(px);
            }
        }

        std::for_each(flaggedPixels.begin(), flaggedPixels.end(), [](auto px) { *px = 0; });

        difference = !flaggedPixels.empty();
        iteration  = !iteration;
    }
}

} // namespace fitting
