"""HTML template constants and CSS for clinical report rendering.

All templates use Python f-string or .format() interpolation.
No external template engines (Jinja2, Mako, etc.) are used.
No em dashes (U+2014) appear in any static text.
"""

# Section marker IDs used for ordering verification in tests.
SECTION_ID_HEADER = "section-header"
SECTION_ID_EXECUTIVE_SUMMARY = "section-executive-summary"
SECTION_ID_FINDINGS_TABLE = "section-findings-table"
SECTION_ID_EVIDENCE_CARDS = "section-evidence-cards"
SECTION_ID_LIMITATIONS = "section-limitations"
SECTION_ID_METHODOLOGY = "section-methodology"
SECTION_ID_SIGN_OFF = "section-sign-off"

# Inlined CSS for the clinical report.
REPORT_CSS = """\
body {
    font-family: "Segoe UI", Arial, sans-serif;
    font-size: 11pt;
    line-height: 1.5;
    color: #1a1a1a;
    margin: 0;
    padding: 2rem;
    max-width: 210mm;
    margin-left: auto;
    margin-right: auto;
}
h1 {
    font-size: 18pt;
    color: #003366;
    border-bottom: 2px solid #003366;
    padding-bottom: 0.3rem;
    margin-top: 1.5rem;
}
h2 {
    font-size: 14pt;
    color: #003366;
    margin-top: 1.5rem;
    margin-bottom: 0.5rem;
}
h3 {
    font-size: 12pt;
    color: #333333;
    margin-top: 1rem;
    margin-bottom: 0.4rem;
}
table {
    width: 100%;
    border-collapse: collapse;
    margin: 1rem 0;
    font-size: 10pt;
}
th {
    background-color: #003366;
    color: #ffffff;
    padding: 0.5rem;
    text-align: left;
    font-weight: 600;
}
td {
    padding: 0.4rem 0.5rem;
    border-bottom: 1px solid #dddddd;
}
tr:nth-child(even) {
    background-color: #f8f9fa;
}
.evidence-card {
    border: 1px solid #cccccc;
    border-radius: 4px;
    padding: 1rem;
    margin: 1rem 0;
    page-break-inside: avoid;
}
.evidence-card h3 {
    margin-top: 0;
}
.summary-box {
    background-color: #f0f4f8;
    border-left: 4px solid #003366;
    padding: 1rem;
    margin: 1rem 0;
}
.limitations-list {
    background-color: #fff8e1;
    border-left: 4px solid #f9a825;
    padding: 1rem;
    margin: 1rem 0;
}
.sign-off-block {
    margin-top: 2rem;
    border-top: 1px solid #cccccc;
    padding-top: 1rem;
}
.sign-off-field {
    margin: 0.5rem 0;
    padding-bottom: 0.3rem;
    border-bottom: 1px dotted #999999;
    min-height: 1.5rem;
}
.metadata {
    color: #666666;
    font-size: 9pt;
}
@media print {
    body {
        padding: 0;
    }
    .evidence-card {
        page-break-inside: avoid;
    }
}
"""

# HTML document skeleton with CSS inlined.
HTML_SKELETON = """\
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clinical Variant Report: {patient_id}</title>
    <style>
{css}
    </style>
</head>
<body>
{body}
</body>
</html>
"""

# Header section template.
HEADER_TEMPLATE = """\
<section id="{section_id}">
    <h1>Clinical Variant Report</h1>
    <div class="metadata">
        <p><strong>Patient ID:</strong> {patient_id}</p>
        <p><strong>Gene Panel:</strong> {panel_name}</p>
        <p><strong>Analysis Date:</strong> {analysis_date}</p>
        <p><strong>Pipeline Version:</strong> {pipeline_version}</p>
    </div>
</section>
"""

# Executive summary section template.
EXECUTIVE_SUMMARY_TEMPLATE = """\
<section id="{section_id}">
    <h2>Executive Summary</h2>
    <div class="summary-box">
        <p><strong>Total variants analyzed:</strong> \
{total_variants_analyzed}</p>
        <p><strong>Variants passed filters:</strong> \
{variants_passed_filters}</p>
        <p><strong>Pathogenic:</strong> {pathogenic_count}</p>
        <p><strong>Likely Pathogenic:</strong> \
{likely_pathogenic_count}</p>
        <p><strong>Variants of Uncertain Significance:</strong> \
{vus_count}</p>
    </div>
</section>
"""

# Findings table section template (header and wrapper).
FINDINGS_TABLE_HEADER = """\
<section id="{section_id}">
    <h2>Findings Table</h2>
    <table>
        <thead>
            <tr>
                <th>Gene</th>
                <th>Consequence</th>
                <th>Classification</th>
                <th>Composite Rank</th>
                <th>Location</th>
            </tr>
        </thead>
        <tbody>
"""

FINDINGS_TABLE_ROW = """\
            <tr>
                <td>{gene_name}</td>
                <td>{consequence}</td>
                <td>{classification}</td>
                <td>{composite_rank}</td>
                <td>{chromosome}:{position}</td>
            </tr>
"""

FINDINGS_TABLE_FOOTER = """\
        </tbody>
    </table>
</section>
"""

# Evidence cards section wrapper.
EVIDENCE_CARDS_HEADER = """\
<section id="{section_id}">
    <h2>Evidence Cards</h2>
"""

EVIDENCE_CARD_TEMPLATE = """\
    <div class="evidence-card">
        <h3>{gene_name}: {consequence}</h3>
        {details}
        <p>{narrative}</p>
    </div>
"""

EVIDENCE_CARDS_FOOTER = """\
</section>
"""

# Evidence card detail line templates.
EVIDENCE_CARD_AF_LINE = """\
        <p><strong>Allele Frequency:</strong> \
{allele_frequency_formatted}</p>
"""

EVIDENCE_CARD_SCORES_LINE = """\
        <p><strong>Predictor Scores:</strong> {scores}</p>
"""

EVIDENCE_CARD_CLINVAR_LINE = """\
        <p><strong>ClinVar:</strong> {clinvar_assertion}</p>
"""

EVIDENCE_CARD_INHERITANCE_LINE = """\
        <p><strong>Inheritance:</strong> {inheritance_pattern}</p>
"""

EVIDENCE_CARD_TAGS_LINE = """\
        <p><strong>ACMG Criteria:</strong> {tags}</p>
"""

# Limitations section template.
LIMITATIONS_TEMPLATE_HEADER = """\
<section id="{section_id}">
    <h2>Limitations</h2>
    <div class="limitations-list">
        <ul>
"""

LIMITATIONS_ITEM = """\
            <li>{limitation}</li>
"""

LIMITATIONS_TEMPLATE_FOOTER = """\
        </ul>
    </div>
</section>
"""

LIMITATIONS_NONE_TEMPLATE = """\
<section id="{section_id}">
    <h2>Limitations</h2>
    <p>No data source limitations were encountered during this analysis.</p>
</section>
"""

# Methodology section template.
METHODOLOGY_TEMPLATE = """\
<section id="{section_id}">
    <h2>Methodology</h2>
    <p><strong>Pipeline Version:</strong> {pipeline_version}</p>
    <p><strong>Analysis Timestamp:</strong> {analysis_timestamp}</p>
    <h3>Reference Files</h3>
    <table>
        <thead>
            <tr>
                <th>File</th>
                <th>SHA-256 Checksum</th>
            </tr>
        </thead>
        <tbody>
{reference_rows}
        </tbody>
    </table>
    <h3>Classification Parameters</h3>
    <table>
        <thead>
            <tr>
                <th>Parameter</th>
                <th>Value</th>
            </tr>
        </thead>
        <tbody>
{parameter_rows}
        </tbody>
    </table>
</section>
"""

METHODOLOGY_REF_ROW = """\
            <tr>
                <td>{path}</td>
                <td><code>{checksum}</code></td>
            </tr>
"""

METHODOLOGY_PARAM_ROW = """\
            <tr>
                <td>{param_name}</td>
                <td>{param_value}</td>
            </tr>
"""

# Sign-off section template.
SIGN_OFF_TEMPLATE = """\
<section id="{section_id}">
    <h2>Sign-off</h2>
    <div class="sign-off-block">
        <div class="sign-off-field">
            <strong>Reviewer:</strong> {reviewer_name}
        </div>
        <div class="sign-off-field">
            <strong>Date:</strong> {review_date}
        </div>
        <div class="sign-off-field">
            <strong>Digital Signature:</strong> {digital_signature}
        </div>
    </div>
</section>
"""
