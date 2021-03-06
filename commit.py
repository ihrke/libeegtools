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

print "Updating configure.in ...",
f = open("configure.in", 'r');
lines=f.readlines();
f.close();
f = open("configure.in", "w");
for line in lines:
    if "AC_INIT([libeegtools]," in line.strip().replace(" ", ""):
        f.write("AC_INIT([libeegtools], [%s], [ihrke@nld.ds.mpg.de])\n"%(nversion));
    else:
        f.write(line);
f.close();
print "Done"

#print "Touching ./doc/mainpage.doc ..."
#status=os.system("touch ./doc/mainpage.doc");
#print "Done (%i)"%status

print "Hit <Enter> to 'cvs commit'";
cmd = "cvs commit"
sys.stdin.readline();
print cmd;
os.system(cmd);


print "Uploading API documentation";
cmd = "sh upload_documentation.sh";
print cmd;
os.system(cmd);
