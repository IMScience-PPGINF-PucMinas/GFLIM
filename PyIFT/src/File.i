
%newobject iftFilename;
%feature("autodoc", "2");
char *iftFilename(const char *pathname, const char *suffix);

%newobject iftDirname;
%feature("autodoc", "2");
char *iftDirname(const char *pathname);

%feature("autodoc", "2");
const char *iftFileExt(const char *pathname);

%newobject iftBasename;
%feature("autodoc", "2");
char *iftBasename(const char *pathname);

