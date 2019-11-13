all: main_server.f
	f90 $$F90FLAGS -o server -s -w main_server.f routs/{eb_v4,iso_routines}.f $$LINK_FNL
