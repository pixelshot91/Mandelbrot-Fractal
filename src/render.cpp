#include "render.hpp"
#include <vector>
#include <cstdint>
#include <cassert>
#include <iostream>
#include <string.h>
#include <thread>
#include <chrono>
#include <cmath>
#include <immintrin.h>
#include <sstream>

//#define CHRONO

using Clock = std::chrono::high_resolution_clock;

struct rgb8_t {
  std::uint8_t r;
  std::uint8_t g;
  std::uint8_t b;
};

rgb8_t heat_lut(float x)
{
  float x0 = 1.f / 4.f;
  float x1 = 2.f / 4.f;
  float x2 = 3.f / 4.f;

  if (x < x0)
  {
    auto g = static_cast<std::uint8_t>(x / x0 * 255);
    return rgb8_t{0, g, 255};
  }
  else if (x < x1)
  {
    auto b = static_cast<std::uint8_t>((x1 - x) / x0 * 255);
    return rgb8_t{0, 255, b};
  }
  else if (x < x2)
  {
    auto r = static_cast<std::uint8_t>((x - x1) / x0 * 255);
    return rgb8_t{r, 255, 0};
  }
  else
  {
    auto b = static_cast<std::uint8_t>((1.f - x) / x0 * 255);
    return rgb8_t{255, b, 0};
  }
}

void render_simd_mt(const unsigned int width,
            const unsigned int height,
            const int n_iterations,
						std::vector<unsigned int>& hist,
						std::vector<int>& nb_iter,
						const unsigned int thread_no, const unsigned int nb_thread,
						const unsigned int aligned_width)
{
	const __m256 xmin = _mm256_set1_ps(-2.5);
	const __m256 ymin = _mm256_set1_ps(-1.0);
	const __m256 xscale = _mm256_set1_ps(3.5/ (width-1));
	const __m256 yscale = _mm256_set1_ps(2.0/ (height-1));
	const __m256 threshold = _mm256_set1_ps(4);
	const __m256 one = _mm256_set1_ps(1);
	const __m256 vec_01234567 = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);

  for (unsigned int y = thread_no; y < (height+1)/2; y+=nb_thread)
	{
		__m256 y_pixel = _mm256_set1_ps(y);
		__m256 y_mb = _mm256_fmadd_ps(y_pixel, yscale, ymin);
		for (unsigned int x = 0; x < width; x+=8)
		{
			__m256 x_pixel = _mm256_add_ps(_mm256_set1_ps(x), vec_01234567);
			__m256 x_mb = _mm256_fmadd_ps(x_pixel, xscale, xmin);
			__m256 zr = _mm256_set1_ps(0);
			__m256 zi = _mm256_set1_ps(0);
			__m256 vec_iter = _mm256_set1_ps(0);
			int it = 0;
			for (it = 0; it < n_iterations; it++) {
				__m256 zr2 = _mm256_mul_ps(zr, zr);
				__m256 zi2 = _mm256_mul_ps(zi, zi);
				__m256 norm2 = _mm256_add_ps(zr2, zi2);
				__m256 mask = _mm256_cmp_ps(norm2, threshold, _CMP_LT_OS);

				// Increase iteration counter of the pixels that hasn't diverged
				vec_iter = _mm256_add_ps(_mm256_and_ps(mask, one), vec_iter);
				if (_mm256_testz_ps(mask, _mm256_set1_ps(-1))) // Every pixel has diverged
					break;

				__m256 zrzi = _mm256_mul_ps(zr, zi);
				zr = _mm256_add_ps(_mm256_sub_ps(zr2, zi2), x_mb);
				zi = _mm256_add_ps(_mm256_add_ps(zrzi, zrzi), y_mb);
			}
			__m256i its_simd = _mm256_cvtps_epi32(vec_iter);

			auto offset = y * aligned_width + x;
			_mm256_storeu_si256((__m256i*)(nb_iter.data() + offset), its_simd);
			
			for (int i = 0; i < 8 && x + i < width; i++) {
				hist[nb_iter[offset + i]] += 2 - (y == height/2);
			}
		}
  }
}

void render(std::byte* buffer,
            const int width,
            const int height,
            const std::ptrdiff_t stride,
            const int n_iterations)
{
  std::vector<unsigned int> hist;
	hist.resize(n_iterations + 1);

	int aligned_width = (width % 8 == 0) ? width : width / 8 * 8 + 8;
	std::vector<int> nb_iter;
	nb_iter.resize(height * aligned_width);

	render_simd_mt(width, height, n_iterations, hist, nb_iter, 0, 1, aligned_width);

	std::vector<rgb8_t> hl;
	hl.resize(n_iterations);
	float total = width * height - hist[n_iterations];	
	hist[0] = 0;
	hl[0] = heat_lut(0);
	int sum = 0;
	for (int i = 1; i <= n_iterations; i++) {
		sum += hist[i];
		hl[i] = heat_lut(sum / total);
	}
  
	for (int y = 0; y < (height+1)/2; y++)
	{
    rgb8_t* lineptr = reinterpret_cast<rgb8_t*>(buffer + stride * y);
    for (int x = 0; x < width; x++) {
			if (nb_iter[y * aligned_width + x] == n_iterations)
				lineptr[x] = {0, 0, 0};
			else
				lineptr[x] = hl[nb_iter[y * aligned_width + x]];
		}
	}
	for (int y = 0; y < height/2; y++) 
		memcpy(buffer + stride * (height - y - 1), buffer + y * stride, stride);
}

void iteration_to_color(std::byte* buffer,
												const std::vector<rgb8_t>& hl,
												const std::vector<int>& nb_iter,
												const int width, const int height, const unsigned int stride,
												const int y_start, const int y_end,
												const int aligned_width)
{
	for (int y = y_start; y < y_end; y++)
	{
    for (int x = 0; x < width; x++)
	  {
			rgb8_t* lineptr = reinterpret_cast<rgb8_t*>(buffer + stride * y);
		  lineptr[x] = hl[nb_iter[y * aligned_width + x]];
		}
		memcpy(buffer + stride * (height - y - 1), buffer + y * stride, stride);
	}
}

void render_mt(std::byte* buffer,
               int width,
               int height,
               std::ptrdiff_t stride,
               int n_iterations)
{
	auto aligned_width = (width % 8 == 0) ? width : width / 8 * 8 + 8;
	std::vector<int> nb_iter;
	nb_iter.resize(height * aligned_width);
	
	unsigned int nb_thread = std::thread::hardware_concurrency();
	std::vector<std::thread> threads_pool;
	
	std::vector<std::vector<unsigned int>> hists;
	hists.resize(nb_thread);

#ifdef CHRONO
	auto start = Clock::now();
#endif
	for (unsigned int n = 0; n < nb_thread; n++) {
		hists[n].resize(n_iterations + 1);
		threads_pool.emplace_back(render_simd_mt, width, height, n_iterations, std::ref(hists[n]), std::ref(nb_iter), n, nb_thread, aligned_width);
	}

	for (auto& t : threads_pool)
		t.join();

#ifdef CHRONO
	std::cout << "Thread total = " << std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start).count() << " ms" << std::endl;
#endif
  
	unsigned int n_max = 0;
	for (unsigned int t = 0; t < hists.size(); t++)
		n_max += hists[t][n_iterations];

	unsigned int total = width * height - n_max;

	std::vector<rgb8_t> hl;
	hl.resize(n_iterations+1);

	unsigned int sum = 0;
	for (int it = 1; it < n_iterations; it++) {
		for (unsigned int t = 0; t < hists.size(); t++)
			sum += hists[t][it];
		hl[it] = heat_lut((float) sum  / total);
	}
	hl[n_iterations] = {0, 0, 0};

#ifdef CHRONO
	std::cout << "HLO = " << std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start).count() << " ms" << std::endl;
#endif
	
	const int nb_line_per_thread = std::ceil(height / 2.0 / nb_thread);
	std::vector<std::thread> threads;
	for (unsigned int n = 0; n < nb_thread; n++) {
		const int y_end = (n == nb_thread - 1) ? height/2 : (n+1)*nb_line_per_thread;
		threads.emplace_back(iteration_to_color, buffer, std::ref(hl), std::ref(nb_iter), width, height, stride, n*nb_line_per_thread, y_end, aligned_width);
	}
	for (auto& t : threads)
		t.join();
	
#ifdef CHRONO
	std::cout << "Total = " << std::chrono::duration_cast<std::chrono::milliseconds>(Clock::now() - start).count() << " ms" << std::endl;
#endif
}
