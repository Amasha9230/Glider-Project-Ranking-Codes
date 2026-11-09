import glob
import os
import re

# 1. DIRECTORY SETUP & FILE OUTPUTS

polar_dir = r"C:\polar_data"  # Folder containing raw XFLR5 polar exports
output_all = os.path.join(polar_dir, "ALL_RANKED_AIRFOILS.csv")  # Destination for full ranking
output_top10 = os.path.join(polar_dir, "TOP_10_WEIGHTED_AIRFOILS.csv")  # Destination for top 10 shortlist


def parse_polar_file(filepath):
  """Reads a single XFLR5 polar file and extracts aerodynamic performance metrics."""
  foil_name = "Unknown"
  re_val = "Unknown"
  records = []

  # Safe file read; skip file if corrupted
  try:
    with open(filepath, "r", errors="ignore") as f:
      lines = f.readlines()
  except Exception:
    return None

  # --- A. Header Extraction ---
  # Search first 15 lines for airfoil name and Reynolds number
  for line in lines[:15]:
    if "Calculated polar for:" in line:
      foil_name = line.replace("Calculated polar for:", "").strip()
    elif "Re =" in line:
      m = re.search(r"Re\s*=\s*([\d\.\s\w]+)", line)
      if m:
        re_val = m.group(1).split()[0].replace(",", "")

  # extract name from filename if missing in header
  if foil_name == "Unknown" or not foil_name:
    base = os.path.basename(filepath).replace(".txt", "").replace(".plr", "")
    foil_name = re.sub(r"_T\d+_Re.*$", "", base)

  # B. Numerical Aerodynamic Data ---
  # Expected columns: Alpha, Cl, Cd, Cdp, Cm
  for line in lines:
    clean = line.strip()
    # Check if row starts with a numeric value
    if clean and re.match(r"^-?\d+\.\d+", clean):
      parts = [float(v) for v in clean.split()]
      if len(parts) >= 5:
        records.append({
            "alpha": parts[0],  # Angle of attack (deg)
            "cl": parts[1],  # Lift coefficient
            "cd": parts[2],  # Drag coefficient
            "cm": parts[4],  # Pitching moment coefficient
        })

  # Discard incomplete or failed simulations
  if len(records) < 4:
    return None

  # --- C. Figure of Merit Tracking Variables ---
  best_ld = -999.0
  best_alpha = 0.0
  best_cl = 0.0
  best_cd = 0.0
  max_cl = -999.0
  max_pf = -999.0
  cm0 = 0.0
  min_alpha_diff = 999.0

  # Scan each angle of attack to evaluate flight metrics
  for r in records:
    a = r["alpha"]
    cl = r["cl"]
    cd = r["cd"]
    cm = r["cm"]

    # Filter for positive physical glide region
    if cd > 0:
      ld = cl / cd
      # Track peak glide ratio (Max distance condition)
      if ld > best_ld:
        best_ld = ld
        best_alpha = a
        best_cl = cl
        best_cd = cd

      # Track loiter / power factor (Min sink speed / max flight time)
      if cl > 0:
        pf = (cl**1.5) / cd
        if pf > max_pf:
          max_pf = pf

    # Track maximum lift coefficient
    if cl > max_cl:
      max_cl = cl

    # Find zero-lift pitching moment (Cm closest to alpha = 0 deg)
    if abs(a - 0.0) < min_alpha_diff:
      min_alpha_diff = abs(a - 0.0)
      cm0 = cm

  # Return single summarized airfoil dictionary
  return {
      "Airfoil": foil_name,
      "Re": re_val,
      "Max_L_D": round(best_ld, 2),
      "Alpha_at_Max_L_D": best_alpha,
      "CL_at_Max_L_D": round(best_cl, 3),
      "CD_at_Max_L_D": round(best_cd, 4),
      "Max_CL": round(max_cl, 3),
      "Max_CL1.5_CD": round(max_pf, 2),
      "Cm0": round(cm0, 3),
      "Valid_Points": len(records),
  }


# ---------------------------------------------------------
# 2. SCAN DIRECTORY & PARSE ALL FILES
# ---------------------------------------------------------
all_files = glob.glob(os.path.join(polar_dir, "*.*"))
results = []

for fp in all_files:
  # Process polar files while skipping output CSVs
  if fp.endswith((".txt", ".plr", ".csv")) and "AIRFOIL" not in fp:
    res = parse_polar_file(fp)
    if res and res["Max_L_D"] > 0:
      results.append(res)

if not results:
  print("No valid polar files found in " + polar_dir)
else:
  # ---------------------------------------------------------
  # 3. DEDUPLICATION
  # ---------------------------------------------------------
  # Keep only the highest Max L/D polar for each unique airfoil name
  unique_foils = {}
  for item in results:
    name = item["Airfoil"]
    if (
        name not in unique_foils
        or item["Max_L_D"] > unique_foils[name]["Max_L_D"]
    ):
      unique_foils[name] = item

  foils_list = list(unique_foils.values())

  # ---------------------------------------------------------
  # 4. NORMALIZATION & COMPOSITE SCORING
  # ---------------------------------------------------------
  # Extract arrays for min-max scaling across all candidates
  all_ld = [f["Max_L_D"] for f in foils_list]
  all_pf = [f["Max_CL1.5_CD"] for f in foils_list]
  all_cm = [abs(f["Cm0"]) for f in foils_list]

  min_ld, max_ld = min(all_ld), max(all_ld)
  min_pf, max_pf = min(all_pf), max(all_pf)
  min_cm, max_cm = min(all_cm), max(all_cm)

  # Scale metrics (0 to 1) and calculate composite utility score
  for f in foils_list:
    n_ld = (f["Max_L_D"] - min_ld) / (
        max_ld - min_ld + 1e-9
    )  # Normalized Glide
    n_pf = (f["Max_CL1.5_CD"] - min_pf) / (
        max_pf - min_pf + 1e-9
    )  # Normalized Loiter
    n_cm = 1.0 - (
        (abs(f["Cm0"]) - min_cm) / (max_cm - min_cm + 1e-9)
    )  # Mild moment reward (inverted penalty)

    # Multi-attribute weights: 50% Distance + 35% Flight Time + 15% Pitch Trim
    score = (0.50 * n_ld + 0.35 * n_pf + 0.15 * n_cm) * 100.0
    f["Score"] = round(score, 1)

  # ---------------------------------------------------------
  # 5. SORTING & LEADERBOARD DISPLAY
  # ---------------------------------------------------------
  # Sort all candidates descending by total score
  ranked = sorted(foils_list, key=lambda x: x["Score"], reverse=True)

  # Print complete leaderboard to terminal
  print("\n" + "=" * 88)
  print(
      f"   COMPLETE AIRFOIL RANKING ({len(ranked)} AIRFOILS): 50% DISTANCE /"
      " 35% FLIGHT TIME"
  )
  print("=" * 88)
  header = f"{'Rank':<5} {'Airfoil':<26} {'Score':<8} {'Max L/D':<9} {'Alpha':<7} {'Cl^1.5/Cd':<11} {'Max Cl':<8} {'Cm0':<6}"
  print(header)
  print("-" * 88)

  for i, f in enumerate(ranked, 1):
    name = f["Airfoil"][:24]
    print(
        f"{i:<5} {name:<26} {f['Score']:<8.1f} {f['Max_L_D']:<9.2f}"
        f" {f['Alpha_at_Max_L_D']:<7.1f} {f['Max_CL1.5_CD']:<11.2f}"
        f" {f['Max_CL']:<8.3f} {f['Cm0']:<6.3f}"
    )
    # Visual marker after Top 10 cutoff
    if i == 10:
      print("-" * 88 + "  [TOP 10 SHORTLIST CUTOFF]")

  print("=" * 88)

  # ---------------------------------------------------------
  # 6. CSV EXPORTS
  # ---------------------------------------------------------
  # Export 1: Complete list of all evaluated airfoils
  with open(output_all, "w") as out:
    out.write(
        "Rank,Airfoil,Score,Max_L_D,Alpha_at_Max_L_D,CL_at_Max_L_D,CD_at_Max_L_D,Max_CL,Max_CL1.5_CD,Cm0,Re\n"
    )
    for i, f in enumerate(ranked, 1):
      out.write(
          f"{i},{f['Airfoil']},{f['Score']},{f['Max_L_D']},{f['Alpha_at_Max_L_D']},{f['CL_at_Max_L_D']},{f['CD_at_Max_L_D']},{f['Max_CL']},{f['Max_CL1.5_CD']},{f['Cm0']},{f['Re']}\n"
      )

  # Export 2: Shortlisted Top 10 airfoils advancing to 3D verification
  with open(output_top10, "w") as out:
    out.write(
        "Rank,Airfoil,Score,Max_L_D,Alpha_at_Max_L_D,CL_at_Max_L_D,CD_at_Max_L_D,Max_CL,Max_CL1.5_CD,Cm0,Re\n"
    )
    for i, f in enumerate(ranked[:10], 1):
      out.write(
          f"{i},{f['Airfoil']},{f['Score']},{f['Max_L_D']},{f['Alpha_at_Max_L_D']},{f['CL_at_Max_L_D']},{f['CD_at_Max_L_D']},{f['Max_CL']},{f['Max_CL1.5_CD']},{f['Cm0']},{f['Re']}\n"
      )

  print(f"\n[OK] Complete list ({len(ranked)} airfoils) saved to: {output_all}")
  print(f"[OK] Top 10 candidates saved to: {output_top10}")
