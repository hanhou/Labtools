CC=g++
CFLAGS=-pg -ggdb -Wall
LFLAGS=-lm -pg
TODAY=`date +%d-%m-%G`
LONGTODAY=`date +%G-%m-%d`
GIT_DESCRIPTION=`git describe --tags`
CLI_VERSION_HEADER=cli_version.h

BUILD=build
SRC=../src

HEADERS= $(addprefix $(SRC)/, core.h data.h errors.h optimizer.h prior.h psychometric.h sigmoid.h bootstrap.h mclist.h special.h mcmc.h rng.h linalg.h getstart.h )
OBJECTS= $(addprefix $(BUILD)/, core.o data.o optimizer.o psychometric.o sigmoid.o bootstrap.o mclist.o special.o mcmc.o rng.o linalg.o getstart.o prior.o)
CLI_H= cli.h cli_utilities.h
CLI_O= $(addprefix $(BUILD)/, cli.o cli_utilities.o)

compile: cli-version psignifit-mapestimate psignifit-bootstrap psignifit-mcmc psignifit-diagnostics

clean:
	-rm -r $(BUILD)

cli-version:
	if git rev-parse &> /dev/null ; then \
		echo "#ifndef CLI_VERSION_H" > $(CLI_VERSION_HEADER) ; \
		echo "#define CLI_VERSION_H" >> $(CLI_VERSION_HEADER) ; \
		echo "#define VERSION \""$(GIT_DESCRIPTION)"\"" >> $(CLI_VERSION_HEADER) ; \
		echo "#endif" >> $(CLI_VERSION_HEADER) ; \
	fi

psignifit-mapestimate: $(OBJECTS) $(HEADERS) $(BUILD)/psignifit-mapestimate.o $(CLI_H) $(CLI_O)
	$(CC) -o psignifit-mapestimate $(OBJECTS) $(BUILD)/psignifit-mapestimate.o $(CLI_O) $(LFLAGS)
psignifit-bootstrap: $(OBJECTS) $(HEADERS) $(BUILD)/psignifit-bootstrap.o $(CLI_H) $(CLI_O)
	$(CC) -o psignifit-bootstrap $(OBJECTS) $(BUILD)/psignifit-bootstrap.o $(CLI_O) $(LFLAGS)
psignifit-mcmc: $(OBJECTS) $(HEADERS) $(BUILD)/psignifit-mcmc.o $(CLI_H) $(CLI_O)
	$(CC) -o psignifit-mcmc $(OBJECTS) $(BUILD)/psignifit-mcmc.o $(CLI_O) $(LFLAGS)
psignifit-diagnostics: $(OBJECTS) $(HEADERS) $(BUILD)/psignifit-diagnostics.o $(CLI_H) $(CLI_O)
	$(CC) -o psignifit-diagnostics $(OBJECTS) $(BUILD)/psignifit-diagnostics.o $(CLI_O) $(LFLAGS)

$(BUILD)/psignifit-mapestimate.o: psignifit-mapestimate.cc $(HEADERS) $(CLI_H)
	$(CC) -c $(CFLAGS) psignifit-mapestimate.cc -o $(BUILD)/psignifit-mapestimate.o
$(BUILD)/psignifit-bootstrap.o: psignifit-bootstrap.cc $(HEADERS) $(CLI_H)
	$(CC) -c $(CFLAGS) psignifit-bootstrap.cc -o $(BUILD)/psignifit-bootstrap.o
$(BUILD)/psignifit-mcmc.o: psignifit-mcmc.cc $(HEADERS) $(CLI_H)
	$(CC) -c $(CFLAGS) psignifit-mcmc.cc -o $(BUILD)/psignifit-mcmc.o
$(BUILD)/psignifit-diagnostics.o: psignifit-diagnostics.cc $(HEADERS) $(CLI_H)
	$(CC) -c $(CFLAGS) psignifit-diagnostics.cc -o $(BUILD)/psignifit-diagnostics.o

$(BUILD)/cli.o: cli.cc cli.h
	$(CC) -c $(CFLAGS) cli.cc -o $(BUILD)/cli.o
$(BUILD)/cli_utilities.o: cli_utilities.cc cli_utilities.h
	$(CC) -c $(CFLAGS) cli_utilities.cc -o $(BUILD)/cli_utilities.o

$(BUILD):
	mkdir $(BUILD)

$(BUILD)/prior.o: $(SRC)/prior.cc $(HEADERS) | $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/prior.cc -o $(BUILD)/prior.o
$(BUILD)/core.o: $(SRC)/core.cc $(HEADERS) | $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/core.cc -o $(BUILD)/core.o
$(BUILD)/data.o: $(SRC)/data.cc $(HEADERS) | $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/data.cc -o $(BUILD)/data.o
$(BUILD)/optimizer.o: $(SRC)/optimizer.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/optimizer.cc -o $(BUILD)/optimizer.o
$(BUILD)/psychometric.o: $(SRC)/psychometric.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/psychometric.cc -o $(BUILD)/psychometric.o
$(BUILD)/sigmoid.o: $(SRC)/sigmoid.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/sigmoid.cc -o $(BUILD)/sigmoid.o
$(BUILD)/mclist.o: $(SRC)/mclist.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/mclist.cc -o $(BUILD)/mclist.o
$(BUILD)/bootstrap.o: $(SRC)/bootstrap.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/bootstrap.cc -o $(BUILD)/bootstrap.o
$(BUILD)/special.o: $(SRC)/special.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/special.cc -o $(BUILD)/special.o
$(BUILD)/mcmc.o: $(SRC)/mcmc.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/mcmc.cc -o $(BUILD)/mcmc.o
$(BUILD)/rng.o: $(SRC)/rng.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/rng.cc -o $(BUILD)/rng.o
$(BUILD)/linalg.o: $(SRC)/linalg.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/linalg.cc -o $(BUILD)/linalg.o
$(BUILD)/getstart.o: $(SRC)/getstart.cc $(HEADERS)| $(BUILD)
	$(CC) -c $(CFLAGS) $(SRC)/getstart.cc -o $(BUILD)/getstart.o

