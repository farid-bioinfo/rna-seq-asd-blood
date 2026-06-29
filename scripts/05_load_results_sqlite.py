import sqlite3
import pandas as pd

# Paths
csv_path = "results/tables/DESeq2_results_family_corrected.csv"
db_path  = "results/deseq2_results.db"

# Load CSV
df = pd.read_csv(csv_path)
print(f"Loaded {len(df)} genes from CSV")

# Connect and write to SQLite
conn = sqlite3.connect(db_path)
df.to_sql("deseq2_results", conn, if_exists="replace", index=False)
print(f"Written to database: {db_path}")

# Test query
query = "SELECT * FROM deseq2_results WHERE padj < 0.05 ORDER BY padj ASC LIMIT 10"
result = pd.read_sql_query(query, conn)
print(f"\nTop hits (padj < 0.05):\n{result}")

conn.close()
