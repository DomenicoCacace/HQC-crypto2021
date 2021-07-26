cmake_minimum_required(VERSION 3.6.0)

project(HQC)

set(BASE_DIR ${CMAKE_CURRENT_LIST_DIR}/src)
enable_language(C ASM)
set(CMAKE_C_STANDARD 99)
set(TARGET_NAME hqc-${SECLVL})

if(${CROSSCOMPILE} STREQUAL "1")

	# Source and header files generated from STM32CubeMX
	set(STM32_DIR ${BASE_DIR}/stm32)
	set(CMSIS_DIR ${STM32_DIR}/driver/cmsis)
	set(HAL_DIR ${STM32_DIR}/driver/hal)

	AUX_SOURCE_DIRECTORY(${HAL_DIR}/library HAL_SOURCES)

	set(SOURCES ${STM32_DIR}/stm32f4xx_hal_msp.c
				${STM32_DIR}/stm32f4xx_it.c
				${STM32_DIR}/syscalls.c
				${STM32_DIR}/system_stm32f4xx.c
				${STM32_DIR}/startup_stm32f401xe.s
				${HAL_SOURCES})

	# Set up the arm-none-eabi- toolchain to cross compile
	set(CMAKE_SYSTEM_PROCESSOR ARM)
	set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
	set(CMAKE_C_COMPILER arm-none-eabi-gcc)
	set(CMAKE_ASM_COMPILER ${CMAKE_C_COMPILER})

	set(CMAKE_SYSROOT ${ARM_TOOLCHAIN_DIR}/../arm-none-eabi)
	set(CMAKE_FIND_ROOT_PATH ${BINUTILS_PATH})
	set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
	set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
	set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

	set(TARGET_NAME ${TARGET_NAME}-ARM)
else()
	set(CMAKE_C_COMPILER gcc)
	set(TARGET_NAME ${TARGET_NAME}-native)
	set(FLAGS "-Wall -pedantic -Wuninitialized -fno-var-tracking-assignments -march=native -O3 -g3 -Wno-unused-function")
endif()

set(SOURCES ${SOURCES}
			${BASE_DIR}/codes/code.c
			${BASE_DIR}/fields/fft.c
			${BASE_DIR}/fields/gf.c
			${BASE_DIR}/fields/gf2x.c
			${BASE_DIR}/hqc/hqc.c
			${BASE_DIR}/hqc/kem.c
			${BASE_DIR}/common/parsing.c
			${BASE_DIR}/codes/reed_muller.c
			${BASE_DIR}/codes/reed_solomon.c
			${BASE_DIR}/common/vector.c
			${BASE_DIR}/lib/fips202.c
			${BASE_DIR}/lib/shake_ds.c
			${BASE_DIR}/lib/shake_prng.c)

set(HEADERS ${HEADERS}
			${BASE_DIR}/codes/code.h
			${BASE_DIR}/fields/fft.h
			${BASE_DIR}/fields/gf.h
			${BASE_DIR}/fields/gf2x.h
			${BASE_DIR}/hqc/hqc.h
			${BASE_DIR}/common/parsing.h
			${BASE_DIR}/codes/reed_muller.h
			${BASE_DIR}/codes/reed_solomon.h
			${BASE_DIR}/common/vector.h
			${BASE_DIR}/lib/domains.h
			${BASE_DIR}/lib/fips202.h
			${BASE_DIR}/lib/shake_ds.h
			${BASE_DIR}/lib/shake_prng.h
			${BASE_DIR}/benchmarking/board_config.h
			${BASE_DIR}/benchmarking/timing_stats.h)

if(${MODE} STREQUAL "FUNCTIONAL")
	add_executable(${TARGET_NAME} ${HEADERS} ${SOURCES} ${BASE_DIR}/benchmarking/functional_test.c)
elseif(${MODE} STREQUAL "TIMING")
	add_executable(${TARGET_NAME} ${HEADERS} ${SOURCES} ${BASE_DIR}/benchmarking/timing_test.c)
else()
	error("PLEASE SPECIFY A TARGET")
endif ()

set(FLAGS "${FLAGS} -DSECURITY_LEVEL=${SECLVL}")

# Set the verbosity level
if(${VERBOSE} STREQUAL "1")
	set(FLAGS "${FLAGS} -DDEBUG")
	set(CMAKE_BUILD_TYPE Debug)
elseif(${VERBOSE} STREQUAL "2")
	set(FLAGS "${FLAGS} -DVERBOSE -DDEBUG")
	set(CMAKE_BUILD_TYPE Debug)
endif()

# Toggle ARM configs
if(${CROSSCOMPILE} STREQUAL "1")
	set(FLAGS "${FLAGS} -DCROSSCOMPILE")
endif()

set_property(TARGET ${TARGET_NAME} APPEND PROPERTY COMPILE_FLAGS ${FLAGS})
target_link_libraries(${TARGET_NAME} m)

if(${CROSSCOMPILE} STREQUAL "1")
	set_target_properties(${TARGET_NAME} PROPERTIES SUFFIX ".elf")

	target_compile_definitions(${TARGET_NAME} PUBLIC
		"__weak=__attribute__((weak))"
		"__packed=__attribute__((__packed__))"
		USE_HAL_DRIVER
		STM32F401xE)

	target_compile_options(${TARGET_NAME} PUBLIC
		-mcpu=cortex-m4
		-mthumb
		-mno-thumb-interwork
		-mfpu=vfp
		-msoft-float
		-Wall
		-ffunction-sections
		-fdata-sections
		-O3)

	set(LINKER_SCRIPT ${STM32_DIR}/STM32F401RETx_FLASH.ld)

	target_link_options(${TARGET_NAME} PUBLIC
		-T ${LINKER_SCRIPT}
		-specs=nosys.specs
		-Wl,--gc-section
		-mcpu=cortex-m4
		-mthumb 
	)

	# Generate assembly listing.
	add_custom_command(
		TARGET ${TARGET_NAME}
		COMMAND "arm-none-eabi-objdump"
		ARGS "-S" "$<TARGET_FILE:${TARGET_NAME}>" ">>" "${TARGET_NAME}.lst")

	# Flash firmware to target.
	find_program(stm32programmer NAMES STM32_Programmer.sh STM32_Programmer_CLI PATHS "/opt/STMicroelectronics/STM32Cube/STM32CubeProgrammer/")
	add_custom_target(flash
		${stm32programmer} -c port=SWD -w "$<TARGET_FILE:${TARGET_NAME}>" -g
		DEPENDS ${TARGET_NAME})
endif()

