# -*- coding: utf-8 -*-
# @Time    : 2021/8/18 10:58
# @Author  : Zhongyi Hua
# @File    : gene_rename.py
# @Usage   :
# @Note    :
# @E-mail  : njbxhzy@hotmail.com

import pandas as pd
from selenium import webdriver
from selenium.webdriver import ActionChains

if __name__ == '__main__':
    driver = webdriver.Chrome()
    driver.get('https://www.plabipd.de/plant_genomes_pa.ep')
    driver.maximize_window()
    driver.implicitly_wait(10)

    elements = driver.find_elements_by_css_selector('body > svg > g > g > a > text[fill="#0000ff"]')
    results = []

    for hover_element in elements:
        _result_dict = {'species': hover_element.text}
        ActionChains(driver).move_to_element(hover_element).perform()
        driver.implicitly_wait(5)
        info_elements = driver.find_elements_by_css_selector('body > div.d3-tip.n > p')
        for info_element in info_elements:
            _result_dict.update(dict(zip(['Date', 'Author', 'Title', 'Journal'], info_element.text.split('\n'))))
        results.append(_result_dict)
    pd.DataFrame(results).to_csv('plant_genome.tsv', sep='\t', index=False)



