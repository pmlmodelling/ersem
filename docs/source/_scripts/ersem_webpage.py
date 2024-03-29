from bs4 import BeautifulSoup
from urllib import request
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

def generator_web_doc(file_name):
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
    output_string = "\n".join(save_string[:1]).replace(
            "ERSEM", "`{0}  <{1}>`__".format("ERSEM", url), 1)

    with open(file_name, "w") as model_file:
        model_file.write(output_string)
