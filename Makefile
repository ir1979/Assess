# Builds all the projects in the solution...
.PHONY: all_projects
all_projects: NPN 

# Builds project 'NPN'...
.PHONY: NPN
NPN: 
	make --directory="NPN/" --file=NPN.makefile

# Cleans all projects...
.PHONY: clean
clean:
	make --directory="NPN/" --file=NPN.makefile clean

