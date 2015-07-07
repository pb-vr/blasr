SHELL=/bin/bash -e -E

.PHONY: curdir submodules p4togit pullfromgit 

# git p4 operations
curdir = $(shell pwd)

submodules = libcpp 

p4togit: $(submodules)
	@for submodule in $(submodules); do \
		cd $$submodule; echo Rebasing submodule $$submodule; git p4 sync; git p4 rebase; git push -u origin master; cd $(curdir);\
	done 
	@for submodule in $(submodules); do \
		echo Adding submodule $$submodle to blasr; git add $$submodule; \
	done
	@git commit -m "Push latest $(submodules) from p4 to github" 
	@git push -u origin master || echo push failed

# To help users sync all submodules from github to local.
pullfromgit:
	@git pull -u origin master
	@git submodule update --init --recursive

