from numpy.f2py.crackfortran import crackfortran
import glob
from contextlib import redirect_stdout
import io


def create(file_name):
    f = io.StringIO()
    main_header = "ERSEM api"
    heading_boarder = "#" * len(main_header)
    api_str = ".. _api::\n\n{}\n{}\n{}\n\n".format(heading_boarder,
                                                      main_header,
                                                      heading_boarder)
    with redirect_stdout(f):
        for fname in glob.glob("../../src/*.F90"):
            if fname != "../../src/shared.F90":
                interface = [sub["saved_interface"].split("\n")[1:-1]
                             for sub in crackfortran(fname)[0]["body"]
                             if sub["name"] != "unknown_type" and sub["block"] in
                             ["subroutine", "function"]]
                subroutines = [".. code-block:: bash \n\n" + "\n".join(i)
                               for i in interface]
                heading = "{} module with the following " \
                              "subroutines".format(fname.split("/")[-1].replace(".F90", ""))
                boarder = "~" * len(heading)
                api_str += "\n\n{}\n{}\n\n{}".format(heading, boarder, "\n".join(subroutines))

    with open(file_name, "w") as api_file:
        api_file.write(api_str)
