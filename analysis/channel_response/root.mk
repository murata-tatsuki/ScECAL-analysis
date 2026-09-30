# Common compiler and ROOT settings for channel_response analyses.
CXX = g++
CXXFLAGS = -O2 -Wall -Wextra -std=c++17
ROOTFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)
