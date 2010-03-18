#!/bin/sh
GREPOPTS=-E
grep $GREPOPTS -n $1 src/*.c
grep $GREPOPTS -n $1 src/*.h
grep $GREPOPTS -n $1 examples/*.c
grep $GREPOPTS -n $1 matlab/*.c
grep $GREPOPTS -n $1 programs/*.c
grep $GREPOPTS -n $1 test/*.c
grep $GREPOPTS -n $1 test/*.h
