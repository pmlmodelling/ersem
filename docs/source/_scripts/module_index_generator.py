from numpy.f2py.crackfortran import crackfortran
from numpy import where
import glob
from contextlib import redirect_stdout
import io
import re


def create(file_name):
    """
    Function to create ERSEM module index everytime the docs are built

    :param file_name: file name of the api
    :type file_name: str
    """
    main_header = "ERSEM module index"
    heading_boarder = "#" * len(main_header)
    index_str = ".. _api:\n\n{}\n{}\n{}\n\n".format(heading_boarder,
                                                    main_header,
                                                    heading_boarder)

    structure = ["include", "!", "module", "extends"]
    f = io.StringIO()
    # removing standard output from console
    with redirect_stdout(f):
        for fname in glob.glob("../../src/*.F90"):
            if fname != "../../src/shared.F90":
                with open(fname, "r") as f:
                    file_contence = f.readlines()
                module_dict = {"include": [], "!": [], "module": [], "extends": []}
                for line in file_contence:
                    check = [s in line for s in structure]
                    index = where(check)[0]
                    if len(index) > 0:
                        key = structure[index[0]]
                        value = line.split(key)[1][1:].strip("\n")
                        if key == "extends":
                            value = value.split(")")[0]
                        module_dict[key].append(value)
                    if "extends" in line:
                        break
                module_dict["class inheritance"] = module_dict.pop("extends")
                module_dict["description"] = module_dict.pop("!")
                module_str = "Module name: {}\n{}\n\n".format(
                        module_dict["module"][0],
                        (len(module_dict["module"][0]) + 13) * "~")
                module_dict.pop('module', None)
                indent = " " * 4
                for key, item in module_dict.items():
                    if key == "description":
                        list_str = \
                                "\n".join(item).replace("\n", "\n" + indent * 2)
                    else:
                        list_str = f"\n{2*indent}.. code-block:: bash \n\n{indent*4}" + \
                                "\n".join(item).replace("\n", "\n" + indent * 4)
                    module_str += f"{indent}**{key.capitalize()}**:\n{indent*2}{list_str}\n\n"

                # Find additional functions within module
                interface = [sub["saved_interface"].split("\n")[1:-1]
                             for sub in crackfortran(fname)[0]["body"]
                             if sub["block"] in ["function"]]
                if interface != []:
                    subroutines = []
                    for i in interface:
                        joined_str = "\n".join(i)
                        subroutines.append(".. code-block:: bash \n\n" + joined_str)
                        func_str = "{}**Functions**:\n\n{}".format(indent,
                                                           "\n".join(subroutines))
                        module_str += func_str.replace("\n", "\n" + indent * 2) +"\n"
            index_str += module_str

        with open(file_name, "w") as index_file:
            index_file.write(index_str)
