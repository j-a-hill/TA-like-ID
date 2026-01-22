#!/usr/bin/env python3
"""
TA-like-ID Shiny Application - Membrane Protein Analysis
Interactive filtering and data exploration with column visibility controls
Designed for Shinylive deployment
"""

from shiny import App, render, ui, reactive
import pandas as pd
import io
from urllib import request

# Load data from GitHub
def load_data():
    """Load the CSV data from GitHub"""
    url = "https://raw.githubusercontent.com/j-a-hill/TA-like-ID/main/raw_data/membrane_protein_analysis_with_reduced_cc.csv"
    try:
        df = pd.read_csv(url)
        return df
    except Exception as e:
        print(f"Error loading data: {e}")
        return pd.DataFrame()

# Load data at startup
df = load_data()

# Get all column names for visibility toggles
ALL_COLUMNS = list(df.columns) if len(df) > 0 else []

# App UI
app_ui = ui.page_fillable(
    ui.tags.style("""
        .shiny-data-grid table tbody tr {
            height: 40px !important;
            max-height: 40px !important;
        }
        .shiny-data-grid table tbody td {
            white-space: nowrap !important;
            overflow: hidden !important;
            text-overflow: ellipsis !important;
            max-height: 40px !important;
        }
        .column-toggles {
            display: flex;
            flex-wrap: wrap;
            gap: 15px;
            padding: 15px;
            background-color: #f8f9fa;
            border-radius: 5px;
            margin-bottom: 15px;
        }
        .filter-group {
            border-bottom: 1px solid #dee2e6;
            padding-bottom: 15px;
            margin-bottom: 15px;
        }
        .filter-group:last-child {
            border-bottom: none;
        }
    """),
    
    ui.card(
        ui.card_header("Column Visibility"),
        ui.output_ui("column_toggles"),
        class_="mb-3"
    ),
    
    ui.layout_sidebar(
        ui.sidebar(
            ui.h4("Filters"),
            
            # Protein Length
            ui.div(
                {"class": "filter-group"},
                ui.input_slider("length_range", "Protein Length", 
                               min=0, max=int(df["Length"].max()) if len(df) > 0 else 1000, 
                               value=[0, int(df["Length"].max()) if len(df) > 0 else 1000])
            ),
            
            # C-terminal Distance
            ui.div(
                {"class": "filter-group"},
                ui.h6("C-term Distance"),
                ui.input_select("cterm_op", "Operator", choices=["All", "<", ">", "="]),
                ui.input_numeric("cterm_val", "Value", value=0, min=0)
            ),
            
            # Penultimate Distance
            ui.div(
                {"class": "filter-group"},
                ui.h6("Penultimate Distance"),
                ui.input_select("penult_op", "Operator", choices=["All", "<", ">", "="]),
                ui.input_numeric("penult_val", "Value", value=0, min=0)
            ),
            
            # Membrane Domain Count
            ui.div(
                {"class": "filter-group"},
                ui.h6("Membrane Domain Count"),
                ui.input_select("md_op", "Operator", choices=["All", "<", ">", "="]),
                ui.input_numeric("md_val", "Value", value=0, min=0)
            ),
            
            # N-terminal Membrane Domains
            ui.div(
                {"class": "filter-group"},
                ui.h6("N-terminal Membrane Domains"),
                ui.input_select("ntmd_op", "Operator", choices=["All", "<", ">", "="]),
                ui.input_numeric("ntmd_val", "Value", value=0, min=0)
            ),
            
            # SRP Prediction
            ui.div(
                {"class": "filter-group"},
                ui.input_select("Prediction", "SRP Prediction", 
                               choices=["ALL", "SP", "OTHER"] if len(df) > 0 else ["ALL"])
            ),
            
            # CC Terms
            ui.div(
                {"class": "filter-group"},
                ui.input_text("cc_search", "CC Terms (search)", placeholder="e.g., ER, secretory")
            ),
            
            ui.download_button("download_csv", "📥 Download CSV", width="100%"),
            
            width=280
        ),
        
        ui.output_data_frame("data_grid")
    )
)


# Server logic
def server(input, output, session):
    
    @output
    @render.ui
    def column_toggles():
        """Create checkboxes for column visibility"""
        checkboxes = []
        for col in ALL_COLUMNS:
            checkboxes.append(
                ui.input_checkbox(f"col_{col.replace('.', '_')}", col, value=True)
            )
        return ui.div({"class": "column-toggles"}, *checkboxes)
    
    @reactive.Calc
    def visible_columns():
        """Get list of visible columns based on checkboxes"""
        visible = []
        for col in ALL_COLUMNS:
            col_id = f"col_{col.replace('.', '_')}"
            if input[col_id]():
                visible.append(col)
        return visible if visible else ALL_COLUMNS
    
    @reactive.Calc
    def filtered_data():
        """Apply all filters to the data"""
        if len(df) == 0:
            return pd.DataFrame()
        
        result = df.copy()
        
        # Length filter
        length_range = input.length_range()
        if length_range:
            result = result[(result["Length"] >= length_range[0]) & (result["Length"] <= length_range[1])]
        
        # C-term distance filter
        cterm_op = input.cterm_op()
        cterm_val = input.cterm_val()
        if cterm_op != "All" and cterm_op in ["<", ">", "="]:
            if cterm_op == "<":
                result = result[result["cterm_distance"] < cterm_val]
            elif cterm_op == ">":
                result = result[result["cterm_distance"] > cterm_val]
            elif cterm_op == "=":
                result = result[result["cterm_distance"] == cterm_val]
        
        # Penultimate distance filter
        penult_op = input.penult_op()
        penult_val = input.penult_val()
        if penult_op != "All" and penult_op in ["<", ">", "="]:
            if penult_op == "<":
                result = result[result["penultimate_distance"] < penult_val]
            elif penult_op == ">":
                result = result[result["penultimate_distance"] > penult_val]
            elif penult_op == "=":
                result = result[result["penultimate_distance"] == penult_val]
        
        # Membrane domain count filter
        md_op = input.md_op()
        md_val = input.md_val()
        if md_op != "All" and md_op in ["<", ">", "="]:
            if md_op == "<":
                result = result[result["membrane_domain_count"] < md_val]
            elif md_op == ">":
                result = result[result["membrane_domain_count"] > md_val]
            elif md_op == "=":
                result = result[result["membrane_domain_count"] == md_val]
        
        # N-terminal membrane domains filter
        ntmd_op = input.ntmd_op()
        ntmd_val = input.ntmd_val()
        if ntmd_op != "All" and ntmd_op in ["<", ">", "="]:
            if ntmd_op == "<":
                result = result[result["N_term_md"] < ntmd_val]
            elif ntmd_op == ">":
                result = result[result["N_term_md"] > ntmd_val]
            elif ntmd_op == "=":
                result = result[result["N_term_md"] == ntmd_val]
        
        # SRP Prediction filter
        pred = input.Prediction()
        if pred and pred != "ALL":
            result = result[result["Prediction"] == pred]
        
        # CC Terms search filter
        cc_search = input.cc_search()
        if cc_search and cc_search.strip():
            search_term = cc_search.strip().lower()
            result = result[result["Reduced.CC.Terms"].str.lower().str.contains(search_term, na=False)]
        
        return result
    
    @output
    @render.data_frame
    def data_grid():
        """Render the filtered data grid with selected columns"""
        fdf = filtered_data()
        visible_cols = visible_columns()
        
        if len(fdf) == 0:
            return pd.DataFrame({"Status": ["Loading data from GitHub..."]})
        
        # Select only visible columns
        display_df = fdf[visible_cols]
        
        return render.DataGrid(
            display_df,
            filters=True,
            width="100%",
            height="600px"
        )
    
    @render.download(filename="filtered_proteins.csv")
    def download_csv():
        """Download filtered data as CSV"""
        fdf = filtered_data()
        csv_buffer = io.StringIO()
        fdf.to_csv(csv_buffer, index=False)
        csv_buffer.seek(0)
        return csv_buffer.getvalue()

# Create the app
app = App(app_ui, server)
