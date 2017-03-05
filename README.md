# drugbank_xml2csv_python
Python script for converting DrugBank XML to relational CSV files. If loading CSVs to a SQL database, a sample query is:  
`select distinct drugbank_id, drugname, inhibitor, antagonist, agonist`  
`from drugbank05_drug2target`  
`join drugbank05_partner_protein using (partner_id)`  
`join drugbank05_drugs using (drugbank_id)`  
`where gene_name = 'KCNH2';`

Tested on DrugBank 5 (3.4.17)

Dependencies: `lxml` 

Outputs:  
- `drugbank05_drugs.csv`
- `drugbank05_drug2target.csv`
- `drugbank05_drug2target_human.csv`
- `drugbank05_drug2enzyme.csv`
- `drugbank05_drug2enzyme_human.csv`
- `drugbank05_drug2transporter.csv`
- `drugbank05_drug2transporter_human.csv`
- `drugbank05_partner_protein.csv`
- `drugbank05_partner_protein_human.csv`
