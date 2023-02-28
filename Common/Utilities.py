import ROOT

class WorkingPointsTauVSmu:
    VLoose = 1
    Loose = 2
    Medium = 3
    Tight = 4

class WorkingPointsTauVSjet:
   VVVLoose = 1
   VVLoose = 2
   VLoose = 3
   Loose = 4
   Medium = 5
   Tight = 6
   VTight = 7
   VVTight = 8

class WorkingPointsTauVSe:
    VVVLoose = 1
    VVLoose = 2
    VLoose = 3
    Loose = 4
    Medium = 5
    Tight = 6
    VTight = 7
    VVTight = 8

class WorkingPointsBoostedTauVSjet:
   VVLoose = 1
   VLoose = 2
   Loose = 3
   Medium = 4
   Tight = 5
   VTight = 6
   VVTight = 7

def ListToVector(list, type="string"):
	vec = ROOT.std.vector(type)()
	for item in list:
		vec.push_back(item)
	return vec

def generate_enum_class(cls):
    enum_string = "enum class {} : int {{\n".format(cls.__name__)
    for name, value in cls.__dict__.items():
        if not name.startswith("__"):
            enum_string += "    {} = {},\n".format(name, value)
    enum_string += "};"
    return enum_string

