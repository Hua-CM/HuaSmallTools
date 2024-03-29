# -*- coding: utf-8 -*-
# @Time    : 2021/8/18 10:58
# @Author  : Zhongyi Hua
# @File    : gene_rename_EVM.py
# @Usage   :
# @Note    :
# @E-mail  : njbxhzy@hotmail.com

import pandas as pd
from selenium import webdriver
from selenium.webdriver import ActionChains

if __name__ == '__main__':
    def per_page(_url):
        driver = webdriver.Chrome()
        driver.get(_url)
        driver.maximize_window()
        driver.implicitly_wait(10)

        elements = driver.find_elements("css selector", 'body > svg > g > g > a > text[fill="#0000ff"]')
        _results = []

        for hover_element in elements:
            _result_dict = {'species': hover_element.text}
            ActionChains(driver).move_to_element(hover_element).perform()
            driver.implicitly_wait(5)
            info_elements = driver.find_elements("css selector", 'body > div.d3-tip.n > p')
            for info_element in info_elements:
                _result_dict.update(dict(zip(['Date', 'Author', 'Title', 'Journal'], info_element.text.split('\n'))))
                _results.append(_result_dict)
        return _results

    results = []
    for _url in ['https://www.plabipd.de/plant_genomes_pn.ep', 'https://www.plabipd.de/plant_genomes_pa.ep']:
        results += per_page(_url)
    pd.DataFrame(results).to_csv('plant_genome.tsv', sep='\t', index=False)



