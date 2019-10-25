import requests
import re
headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:65.0) Gecko/20100101 Firefox/65.0"}
file_name = r"E:\OtherData\pueraria\SRX.txt"
out_file = r"E:\OtherData\pueraria\SRR.txt"


def get_html(url):
    resp = requests.get(url, headers=headers)
    return resp.text


with open(file_name, "r") as num_list:
    a = re.compile(r"[D,E,S]RR(\d{6,8})")
    for line in num_list:
        tmp_text = get_html("https://www.ncbi.nlm.nih.gov/sra/?term=" + line.strip())
        b = re.search(a, tmp_text).group()
        c = "/".join(["anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra", b[0:3], b[0:6], b, b+".sra"])
        with open(out_file, "a") as f:
            f.writelines(line.strip()+"\t"+b+"\t"+c+"\n")

