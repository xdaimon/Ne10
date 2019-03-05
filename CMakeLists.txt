cmake_minimum_required(VERSION 2.6)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpu=neon -Wno-macro-redefined -Wno-deprecated -fintegrated-as -fno-strict-aliasing")

if(${NE10_ANDROID_TARGET_ARCH} STREQUAL "armv7")
    set(NE10_TARGET_ARCH "armv7")
else()
    set(NE10_TARGET_ARCH "aarch64")
endif()

# Adding cflags for armv7. Aarch64 does not need such flags.
if(${NE10_TARGET_ARCH} STREQUAL "armv7")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -marm -march=armv7-a -mfloat-abi=softfp -Wl,--no-warn-mismatch")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_C_FLAGS}")

add_definitions(-DNE10_UNROLL_LEVEL=1)
include_directories (
    inc
    common
)
add_library(NE10 STATIC
    modules/dsp/NE10_fft_float32.neonintrinsic.c
    modules/dsp/NE10_rfft_float32.neonintrinsic.c
    modules/dsp/NE10_fft_generic_float32.neonintrinsic.cpp
    modules/dsp/NE10_fft.c
    modules/dsp/NE10_rfft_float32.c
)

if (RELEASE_BUILD)
    target_compile_options(NE10 PUBLIC "$<$<CONFIG:RELEASE>:-Ofast>")
    target_compile_options(NE10 PUBLIC "$<$<CONFIG:DEBUG>:-Ofast>")
endif()
target_compile_options(NE10 PUBLIC "$<$<CONFIG:RELEASE>:-mfpu=neon>")
target_compile_options(NE10 PUBLIC "$<$<CONFIG:DEBUG>:-mfpu=neon>")

