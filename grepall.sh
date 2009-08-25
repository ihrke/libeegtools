#!/bin/sh

grep -n $1 src/*.c
grep -n $1 src/*.h
grep -n $1 examples/*.c
grep -n $1 matlab/*.c
grep -n $1 programs/*.c
