from bs4 import BeautifulSoup
from urllib import request
import pdfplumber
import os
import pandas


def generator_web_doc(file_name):
    main_header = "The ERSEM Model"
    heading_boarder = "#" * len(main_header)
    output_string = ".. _model:\n\n{}\n{}\n{}\n\n".format(heading_boarder,
                                                          main_header,
                                                          heading_boarder)

    url = "https://www.pml.ac.uk/Modelling_at_PML/Models/ERSEM"
    html = request.urlopen(url).read()

    soup = BeautifulSoup(html, 'html.parser')
    base_html = \
        "https://www.pml.ac.uk/Modelling/Models/Physical_models_and_couplers#"
    save_string = []
    for string in soup.find_all("p"):
        string = string.get_text().replace("\xa0", " ")
        for software in ["GOTM", "NEMO", "FVCOM", "FABM"]:
            rst_link = "`{0} <{1}{0}>`__".format(software, base_html)
            string = string.replace(software, rst_link)
        save_string.append(string)
    output_string += "\n".join(save_string[:1]).replace(
            "ERSEM", "`{0}  <{1}>`__".format("ERSEM", url), 1)

    with open(file_name, "w") as model_file:
        model_file.write(output_string)


def generator_history(file_name):
    url = \
        "https://www.pml.ac.uk/getmedia/2bf146b7-17d6-41d3-9047-fd17b773f749/ERSEM_history_2.pdf"
    temp_file_name = "temp.pdf"
    response = request.urlretrieve(url, temp_file_name)
    main_header = "The History of ERSEM"
    heading_boarder = "#" * len(main_header)
    output_string = ".. _history:\n\n{}\n{}\n{}\n\n".format(heading_boarder,
                                                            main_header,
                                                            heading_boarder)
    output_page = []
    strings_to_remove = [main_header]
    with pdfplumber.open(temp_file_name) as pdf:
        for page in pdf.pages:
            page_txt = page.extract_text()
            for string in strings_to_remove:
                page_txt = page_txt.replace(string, "")
            output_page.append(page_txt)
    figure = "\n\n.. figure:: ../../images/ERSEM.png\n   :alt: ERSEM diagram\n    :width: 100.0%"
    output_string += "\n".join(output_page) + figure
    with open(file_name, "w") as history_file:
        history_file.write(output_string)
    os.remove(temp_file_name)

