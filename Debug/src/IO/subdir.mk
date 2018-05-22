################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/IO/FastaParser.cpp \
../src/IO/FastqParser.cpp 

OBJS += \
./src/IO/FastaParser.o \
./src/IO/FastqParser.o 

CPP_DEPS += \
./src/IO/FastaParser.d \
./src/IO/FastqParser.d 


# Each subdirectory must supply rules for building sources it contributes
src/IO/%.o: ../src/IO/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


