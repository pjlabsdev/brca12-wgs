//! Per cancer type: which of BRCA1/BRCA2 is relevant and what the result of this analysis means
//! (port of the 1.1 cancer-risks module). Associations and risk figures with their sources come
//! from assets/cancer_associations.json (cohort studies with PubMed IDs). Research use only: the
//! statements describe what was found in the data, not a personal risk.
use serde_json::{Value, json};

const TABLE: &str = include_str!("../assets/cancer_associations.json");

fn interpretation(gene: &str, association: &str, status: &str, risk: &Value) -> Value {
    let r = |l: &str| risk[l].as_str().unwrap_or("").to_string();
    if association == "not_established" {
        return json!({"en": format!("{gene} is not an established risk gene for this cancer type; the result does not change that."),
                      "de": format!("{gene} ist kein gesichertes Risikogen für diese Krebsart; das Ergebnis ändert daran nichts.")});
    }
    match status {
        "pathogenic" => json!({
            "en": format!("Pathogenic variant in {gene} found: risk for carriers {}. Confirmation in an accredited laboratory and genetic counselling are needed.", r("en")),
            "de": format!("Pathogene Variante in {gene} gefunden: Risiko für Trägerinnen und Träger {}. Bestätigung in einem akkreditierten Labor und humangenetische Beratung sind nötig.", r("de"))}),
        "candidate" => json!({
            "en": format!("Finding in {gene} that needs confirmation; if confirmed as pathogenic, the risk for carriers is {}.", r("en")),
            "de": format!("Befund in {gene}, der bestätigt werden muss; falls er sich als pathogen bestätigt, liegt das Risiko für Trägerinnen und Träger bei {}.", r("de"))}),
        "not_assessable" => json!({
            "en": format!("{gene} not assessable: the exons are covered too little for a statement."),
            "de": format!("{gene} nicht beurteilbar: die Exons sind für eine Aussage zu gering abgedeckt.")}),
        _ => {
            let (qe, qd) = if association == "possible" { (" (association not certain)", " (Zusammenhang nicht gesichert)") } else { ("", "") };
            json!({"en": format!("No pathogenic {gene} variant found: no {gene}-related increased risk detected{qe}."),
                   "de": format!("Keine pathogene {gene}-Variante gefunden: kein {gene}-bedingt erhöhtes Risiko erkennbar{qd}.")})
        }
    }
}

/// cancer_risks block of the result JSON; `status` gives the gene result (pathogenic, candidate,
/// none, not_assessable) of BRCA1 and BRCA2
pub fn risks(status: &[(&str, &str)]) -> Value {
    let t: Value = serde_json::from_str(TABLE).expect("bundled cancer table");
    let st = |g: &str| status.iter().find(|x| x.0 == g).map(|x| x.1).unwrap_or("none");
    let cancers: Vec<Value> = t["cancers"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| {
            let genes: Vec<Value> = ["BRCA1", "BRCA2"]
                .iter()
                .map(|g| {
                    let info = &c["genes"][*g];
                    let assoc = info["association"].as_str().unwrap_or("");
                    json!({"gene": g, "association": assoc, "association_label": t["associationLabels"][assoc],
                           "risk_for_carriers": info["risk"], "sources": info["sources"], "status": st(g),
                           "interpretation": interpretation(g, assoc, st(g), &info["risk"])})
                })
                .collect();
            let mut e = json!({"id": c["id"], "name": c["name"], "genes": genes, "other_genes_not_analysed": c["otherGenes"]});
            if !c["note"].is_null() {
                e["note"] = c["note"].clone();
                e["note_sources"] = c["noteSources"].clone();
            }
            e
        })
        .collect();
    json!({
        "version": t["version"], "cancers": cancers, "other_genes_note": t["otherGenesNote"], "sources": t["sources"],
        "disclaimer": {
            "en": "Research use only. The statements describe what this analysis found in BRCA1 and BRCA2. They are not a personal risk estimate and do not cover other genes, regions with low coverage or variant types the method cannot see.",
            "de": "Nur für Forschung. Die Aussagen beschreiben, was diese Analyse in BRCA1 und BRCA2 gefunden hat. Sie sind keine persönliche Risikoeinschätzung und gelten nicht für andere Gene, Bereiche mit geringer Abdeckung oder Variantentypen, die die Methode nicht erfasst.",
        }
    })
}
