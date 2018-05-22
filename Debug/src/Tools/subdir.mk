################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Tools/Filter_length.cpp \
../src/Tools/Nucleotide_distribution.cpp \
../src/Tools/Summarize_length.cpp 

OBJS += \
./src/Tools/Filter_length.o \
./src/Tools/Nucleotide_distribution.o \
./src/Tools/Summarize_length.o 

CPP_DEPS += \
./src/Tools/Filter_length.d \
./src/Tools/Nucleotide_distribution.d \
./src/Tools/Summarize_length.d 


# Each subdirectory must supply rules for building sources it contributes
src/Tools/%.o: ../src/Tools/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


