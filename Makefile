CXX = g++
CXXFLAGS = -g -std=c++20 -Wall -Wextra
LDFLAGS = -lsfml-graphics -lsfml-window -lsfml-system

SRC = main.cpp
OBJDIR = build
OBJ = $(SRC:%.cpp=$(OBJDIR)/%.o)
DEP = $(OBJ:.o=.d)
TARGET = app

# Default target
all: $(TARGET)

# Create build dir if needed
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Rule to compile .cpp → build/.o with dependency generation
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Link all .o into final executable
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@ $(LDFLAGS)

# Include dependency files if they exist
-include $(DEP)

# Cleanup
clean:
	rm -rf $(OBJDIR) $(TARGET)
