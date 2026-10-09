########################################################################
#                                                                      #
#             Makefile for Scripts                                     #
#                                                                      #
#             03-04-2007                                               #
#                                                                      #
########################################################################
#
# Default target:

all: Cinc_field Cscale_field

# Compiler

CC = cc

# Scripts written in C:

Cinc_field: Cinc_field.c
	$(CC) -o Cinc_field Cinc_field.c -lm

Cscale_field: Cscale_field.c
	$(CC) -o Cscale_field Cscale_field.c -lm

# Special targets:

new:
	touch Cinc_field.c Cscale_field.c
	make
