######################################################################
#
#   See COPYING file distributed along with the psignifit package for
#   the copyright and license terms
#
######################################################################

# Makefile for targets in the src subdirectory

CC=g++
CFLAGS=-pg -ggdb -Wall
LFLAGS=-lm -pg

BUILD=build
HEADERS=core.h data.h errors.h optimizer.h prior.h psychometric.h sigmoid.h bootstrap.h mclist.h special.h mcmc.h rng.h linalg.h getstart.h integrate.h
OBJECTS= $(addprefix $(BUILD)/, core.o data.o optimizer.o psychometric.o sigmoid.o bootstrap.o mclist.o special.o mcmc.o rng.o linalg.o getstart.o prior.o integrate.o)
TESTS=tests_all

tests_all: tests_all.cc $(HEADERS) $(OBJECTS) testing.h
	$(CC) -c $(CFLAGS) tests_all.cc -o $(BUILD)/tests_all.o
	$(CC) $(OBJECTS) $(BUILD)/tests_all.o $(LFLAGS) -o tests_all

test: tests_all
	./tests_all

play: $(OBJECTS) play.cc $(HEADERS)
	$(CC) -c $(CFLAGS) play.cc -o $(BUILD)/play.o
	$(CC) $(OBJECTS) $(BUILD)/play.o $(LFLAGS) -o play

$(BUILD):
	mkdir $(BUILD)

$(BUILD)/prior.o: prior.cc $(HEADERS) | $(BUILD)
	$(CC)  -c $(CFLAGS) prior.cc -o $(BUILD)/prior.o
$(BUILD)/core.o: core.cc $(HEADERS) | $(BUILD)
	$(CC) -c $(CFLAGS) core.cc -o $(BUILD)/core.o
$(BUILD)/data.o: data.cc $(HEADERS) | $(BUILD)
	$(CC) -c $(CFLAGS) data.cc -o $(BUILD)/data.o
$(BUILD)/optimizer.o: optimizer.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) optimizer.cc -o $(BUILD)/optimizer.o
$(BUILD)/psychometric.o: psychometric.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) psychometric.cc -o $(BUILD)/psychometric.o
$(BUILD)/sigmoid.o: sigmoid.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) sigmoid.cc -o $(BUILD)/sigmoid.o
$(BUILD)/mclist.o: mclist.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) mclist.cc -o $(BUILD)/mclist.o
$(BUILD)/bootstrap.o: bootstrap.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) bootstrap.cc -o $(BUILD)/bootstrap.o
$(BUILD)/special.o: special.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) special.cc -o $(BUILD)/special.o
$(BUILD)/mcmc.o: mcmc.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) mcmc.cc -o $(BUILD)/mcmc.o
$(BUILD)/rng.o: rng.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) rng.cc -o $(BUILD)/rng.o
$(BUILD)/linalg.o: linalg.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) linalg.cc -o $(BUILD)/linalg.o
$(BUILD)/getstart.o: getstart.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) getstart.cc -o $(BUILD)/getstart.o
$(BUILD)/integrate.o: integrate.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) integrate.cc -o $(BUILD)/integrate.o

clean:
	-rm -rf $(BUILD)
	-rm tests_all
