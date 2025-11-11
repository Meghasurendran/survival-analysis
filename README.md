## 📊 Comparison of the Cox PH Model and Lasso-Based Cox Regression Model: An Application of Head and Neck Cancer

This project compares the **Cox Proportional Hazards (PH) model** and the **Lasso-based Cox regression model** to identify key prognostic factors influencing the survival of **Head and Neck Cancer (HNC)** patients.

### 🏥 Dataset  
The analysis was based on data from **972 patients** treated at **Malabar Cancer Centre, Thalassery**, with a **median follow-up time of 72 months** and a **mean age of 58.84 years**.  
- **Deaths:** 40.43%  
- **Right-censored cases:** 59.57%  
- **Variables included:** age, sex, tumour site, histology, tumour grade, composite stage, tumour size, nodal status, metastasis status, treatment type, and lifestyle factors such as smoking, alcohol consumption, and tobacco chewing.

### 🔍 Methodology  
- **Kaplan-Meier method** was used to estimate survival probabilities.  
- **Log-Rank test** identified significant survival-related factors.  
- **Cox PH** and **Lasso-based Cox regression models** were compared based on their predictive accuracy.  
- Regularization in the Lasso model helped in variable selection and model optimization.

### 📈 Key Findings  
- Both models showed **similar predictive accuracy**:  
  - **C-index:** Cox = 0.676 | Lasso = 0.673  
  - **AUC:** Cox = 0.72 | Lasso = 0.717  
- Significant predictors (Cox PH): tumour site, grade, size, nodal status, treatment type, and smoking habit.  
- Significant predictors (Lasso Cox): age category, primary site, grade, tumour size, nodal status, metastasis, stage, treatment type, and lifestyle factors.  
- Both models demonstrated strong potential for **risk prediction** in head and neck cancer survival analysis.

### 🧰 Tools Used  
- **R** – Survival analysis, model comparison, visualization  
- **SPSS** – Data preprocessing and statistical validation
