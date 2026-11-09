#pragma once
#include <TASImage.h>
#include <TImageDump.h>
#include <TVirtualPad.h>
#include <algorithm>
#include <memory>
#include <cstring>
#include <chrono>
#include <stdexcept>

// ROOT's circle painter scans a whole image for each marker. Paint exactly
// the same circle over a small copy of its background and copy it back.
class BoundedCircleImage : public TASImage {
    std::unique_ptr<TASImage> tile;
public:
    void DrawCircle(Int_t x, Int_t y, Int_t radius, const char* color, Int_t thick) override {
        if (radius < 0 || radius > 100 || thick < -100 || thick > 100) {
            TASImage::DrawCircle(x, y, radius, color, thick);
            return;
        }
        const int reach = radius + std::abs(thick) + 4;
        const int left = std::max(0, x - reach), top = std::max(0, y - reach);
        const int right = std::min(int(GetWidth()), x + reach + 1);
        const int bottom = std::min(int(GetHeight()), y + reach + 1);
        if (left >= right || top >= bottom) return;
        const unsigned width = right - left, height = bottom - top;
        if (!tile || tile->GetWidth() != width || tile->GetHeight() != height)
            tile = std::make_unique<TASImage>(width, height);
        auto* pixels = GetArgbArray();
        auto* local = tile->GetArgbArray();
        if (!pixels || !local) throw std::runtime_error("Cannot allocate marker raster");
        for (unsigned row = 0; row < height; ++row)
            std::memcpy(local + row * width, pixels + (top + row) * GetWidth() + left, width * sizeof(UInt_t));
        tile->DrawCircle(x - left, y - top, radius, color, thick);
        local = tile->GetArgbArray();
        for (unsigned row = 0; row < height; ++row)
            std::memcpy(pixels + (top + row) * GetWidth() + left, local + row * width, width * sizeof(UInt_t));
    }
};

class FastImageDump : public TImageDump {
public:
    FastImageDump() {
        fImage = new BoundedCircleImage;
        fType = 114; // ROOT's in-memory image mode: destructor must not write.
        SetBit(BIT(11)); // kPrintingPS, as in TASImage::FromPad.
    }
};

struct PngTiming {double raster, encode;};
inline PngTiming fastPng(TVirtualPad& canvas, const char* filename) {
    const auto start = std::chrono::steady_clock::now();
    struct Restore {TVirtualPS* old = gVirtualPS; ~Restore() {gVirtualPS = old;}} restore;
    FastImageDump painter;
    TVirtualPad::TContext context(&canvas, kFALSE);
    canvas.Paint();
    auto* image = painter.GetImage();
    if (!image || !image->IsValid()) throw std::runtime_error("Cannot rasterize canvas");
    const auto painted = std::chrono::steady_clock::now();
    image->SetImageCompression(10); // ROOT's 0..100 scale: lossless PNG level 1.
    image->WriteImage(filename, TImage::kPng);
    return {std::chrono::duration<double>(painted-start).count(),
        std::chrono::duration<double>(std::chrono::steady_clock::now()-painted).count()};
}
