from numpy.f2py.crackfortran import crackfortran
import glob
from contextlib import redirect_stdout
import io


def create(file_name):
    """
    Function to create ERSEM api everytime the docs are built

    :param file_name: file name of the api
    :type file_name: str
    """
    main_header = "ERSEM api"
    heading_boarder = "#" * len(main_header)
    api_str = ".. _api:\n\n{}\n{}\n{}\n\n".format(heading_boarder,
                                                      main_header,
                                                      heading_boarder)
    f = io.StringIO()
    # removing standard output from console
    with redirect_stdout(f):
        for fname in glob.glob("../../src/*.F90"):
            if fname != "../../src/shared.F90":
                interface = [sub["saved_interface"].split("\n")[1:-1]
                             for sub in crackfortran(fname)[0]["body"]
                             if sub["name"] != "unknown_type" and sub["block"] in
                             ["subroutine", "function"]]
                subroutines = []
                for i in interface:
                    joined_str = "\n".join(i)
                    multi_subroutine = joined_str.find("register_bn")
                    if multi_subroutine == -1:
                        subroutines.append(".. code-block:: bash \n\n" + joined_str)
                    else:
                        for string in [joined_str[:multi_subroutine - 15],
                                       joined_str[multi_subroutine - 15:-30]]:
                            subroutines.append(".. code-block:: bash \n\n" + string)
                heading = "{} module with the following " \
                              "subroutines".format(fname.split("/")[-1].replace(".F90", ""))
                boarder = "~" * len(heading)
                api_str += "\n\n{}\n{}\n\n{}".format(heading, boarder, "\n".join(subroutines))
    api_str = api_str.replace("do_bn", "do")
    with open(file_name, "w") as api_file:
        api_file.write(api_str)
