# Define compiler
CC = gcc

# Define compiler flags
CFLAGS = -Wall -g

# Define libraries to link
LIBS = -lm

# Define the target executable
TARGET = simulation

# Define the source files
SRCS = simulation.c model.c

# Define the object files
OBJS = $(SRCS:.c=.o)

# Default target
all: $(TARGET)

# Rule to build the target
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Rule to build object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up the build files
clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean