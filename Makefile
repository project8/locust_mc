MONARCHLOC := /usr/local
CFLAGS := -Wall -I $(MONARCHLOC)/include -g
#CFLAGS := -Wall -I $(MONARCHLOC)/include -O3
LIBS := -L$(MONARCHLOC)/lib -lMonarchCore -lfftw3f -lfftw3f_threads -lpthread -lm -lprotobuf -lyajl

generate_data: generate_data.cc
	g++ $(CFLAGS) -c generate_data.cc -o generate_data.o
	g++ $(CFLAGS) -o generate_data $(LIBS) generate_data.o

check_egg_digitization: check_egg_digitization.cc
	g++ $(CFLAGS) -c check_egg_digitization.cc -o check_egg_digitization.o
	g++ $(CFLAGS) -o check_egg_digitization $(LIBS) check_egg_digitization.o
