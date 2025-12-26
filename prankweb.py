"""
PrankWeb binding site prediction module
"""

import os
import time
import zipfile
import pandas as pd
import gradio as gr
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from config import current_pdb_info, PRANKWEB_OUTPUT_DIR


def run_prankweb_prediction():
    """Run PrankWeb prediction on the current PDB file."""
    if not current_pdb_info["pdb_id"] or not current_pdb_info["pdb_path"]:
        return (
            gr.update(value="<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ No structure loaded. Please search for a disease first.</div>", visible=True),
            gr.update(value=None, visible=False)
        )
    
    # Show processing message
    yield (
        gr.update(value="<div style='padding: 20px; background: #fff3cd; border-radius: 8px; color: #856404;'>🔮 Processing PrankWeb prediction (this may take several minutes)...</div>", visible=True),
        gr.update(value=None, visible=False)
    )
    
    pdb_path = current_pdb_info["pdb_path"]
    pdb_id = current_pdb_info["pdb_id"]
    output_dir = PRANKWEB_OUTPUT_DIR
    os.makedirs(output_dir, exist_ok=True)
    
    absolute_path = os.path.abspath(pdb_path)
    
    # Setup Chrome driver with download preferences and HEADLESS mode
    chrome_options = webdriver.ChromeOptions()
    chrome_options.add_argument('--headless=new')
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')
    chrome_options.add_argument('--window-size=1920,1080')
    prefs = {
        "download.default_directory": os.path.abspath(output_dir),
        "download.prompt_for_download": False,
    }
    chrome_options.add_experimental_option("prefs", prefs)
    
    try:
        driver = webdriver.Chrome(options=chrome_options)
        
        driver.get("https://prankweb.cz/")
        time.sleep(3)
        
        wait = WebDriverWait(driver, 30)
        custom_structure = wait.until(EC.presence_of_element_located((By.XPATH, "//*[contains(text(), 'Custom structure')]")))
        driver.execute_script("arguments[0].click();", custom_structure)
        time.sleep(1)
        
        file_input = driver.find_element(By.CSS_SELECTOR, "input[type='file']")
        file_input.send_keys(absolute_path)
        time.sleep(2)
        
        submit_btn = wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, "button[type='submit']")))
        driver.execute_script("arguments[0].click();", submit_btn)
        
        wait_long = WebDriverWait(driver, 600)
        info_tab = wait_long.until(EC.presence_of_element_located((By.XPATH, "//*[contains(text(), 'Info')]")))
        
        driver.execute_script("arguments[0].click();", info_tab)
        time.sleep(2)
        
        download_btn = wait_long.until(EC.presence_of_element_located((By.XPATH, "//*[contains(text(), 'Download prediction data')]")))
        driver.execute_script("arguments[0].click();", download_btn)
        
        time.sleep(10)
        driver.quit()
        
        # Find and extract the zip file
        zip_files = [f for f in os.listdir(output_dir) if f.endswith('.zip')]
        if not zip_files:
            return (
                gr.update(value="❌ Download failed - no zip file found", visible=True),
                gr.update(value=None, visible=False)
            )
        
        zip_path = os.path.join(output_dir, zip_files[0])
        extract_path = os.path.join(output_dir, zip_files[0].replace('.zip', ''))
        
        os.makedirs(extract_path, exist_ok=True)
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
        
        # Find the CSV file
        csv_path = os.path.join(extract_path, "structure.pdb_predictions.csv")
        if not os.path.exists(csv_path):
            return (
                gr.update(value="❌ CSV file not found in extracted data", visible=True),
                gr.update(value=None, visible=False)
            )
        
        # Read and filter CSV
        df = pd.read_csv(csv_path)
        columns_to_drop = ['residue_ids', 'surf_atom_ids']
        df = df.drop(columns=[col for col in columns_to_drop if col in df.columns], errors='ignore')
        
        # Store CSV path globally for later use in docking
        current_pdb_info["prankweb_csv"] = csv_path
        
        yield (
            gr.update(value="<div style='padding: 20px; background: #d4edda; border-radius: 8px; color: #155724;'>✅ PrankWeb prediction completed!</div>", visible=True),
            gr.update(value=df, visible=True)
        )
        
    except Exception as e:
        if 'driver' in locals():
            driver.quit()
        yield (
            gr.update(value=f"<div style='padding: 20px; background: #fee; border-radius: 8px; color: #c33;'>❌ Error: {str(e)}</div>", visible=True),
            gr.update(value=None, visible=False)
        )
