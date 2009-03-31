#!/usr/bin/env python
import os,sys;

f = open("VERSION", "r");
version = f.read().strip();
f.close();

(major,minor,revision) = version.split('.');
nversion = "%i.%i.%i"%(int(major),int(minor),int(revision)+1);
print "Update from version '%s' to '%s'"%(version, nversion)

print "Update Doxyfile ...",
f = open("Doxyfile", 'r');
lines=f.readlines();
f.close();
f = open("Doxyfile", "w");
for line in lines:
    if "PROJECT_NUMBER=" in line.strip().replace(" ", ""):
        f.write("PROJECT_NUMBER = %s\n"%(nversion));
    else:
        f.write(line);
f.close();
print "Done"


print "Updating VERSION ...",
f = open("VERSION", "w");
f.write(nversion);
f.close()
print "Done"

#print "Touching ./doc/mainpage.doc ..."
#status=os.system("touch ./doc/mainpage.doc");
#print "Done (%i)"%status

cmd = "cvs commit"
sys.stdin.readline();
print cmd;
os.system(cmd);
