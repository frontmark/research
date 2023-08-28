ifeq ($(OS),Windows_NT)
ppth = $(shell cygpath -a -m -l .)
else
ppth = $(shell pwd -P)
endif

name = gen_img_from_lay
img  = frontmark/$(name):latest

caps = --cap-drop ALL
sopt = --security-opt no-new-privileges
rofs = --read-only --tmpfs /tmp
vol  = --mount type=bind,source='$(ppth)',target=/home/snake
args = $(caps) $(sopt) $(rofs) $(vol)

.PHONY: all build run stop rm up down init reinit rerun
.DEFAULT_GOAL := all
all: build

build:
	docker build --pull -t $(img) .
run:
	docker run -it --rm $(args) --name $(name)1 --hostname $(name)1 $(img)
stop rm:
	-docker container $@ $(name)1

up: build run

down: stop rm

init: up

reinit: down init

rerun: down run
