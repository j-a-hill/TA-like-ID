#!/usr/bin/env python3
"""
TA-like-ID Membrane Protein Analysis - Shiny App
Interactive filtering and visualization

Designed for easy sharing via Gist or direct deployment.
Loads data from GitHub (works with Shinylive browser deployment).
"""

from shiny import App, render, ui, reactive, Inputs, Outputs, Session
from htmltools import Tag
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import io
from pathlib import Path

# Load data - GitHub first, local fallback
DATA_URL = "https://raw.githubusercontent.com/j-a-hill/TA-like-ID/main/raw_data/membrane_protein_analysis_with_reduced_cc.csv"
_local_file = Path(__file__).parent.parent / 'raw_data' / 'membrane_protein_analysis_with_reduced_cc.csv'

try:
    df = pd.read_csv(DATA_URL)
except Exception:
    try:
        df = pd.read_csv(str(_local_file))
    except Exception as e:
        raise Exception(f"Could not load data: {e}")

# App UI
app_ui = ui.page_fluid(
    ui.tags.head(
        ui.tags.style("""
            body {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                min-height: 100vh;
                padding: 20px;
            }
            .container-fluid {
                background: white;
                border-radius: 12px;
                box-shadow: 0 20px 60px rgba(0,0,0,0.3);
                padding: 0;
                overflow: hidden;
            }
            .app-header {
                background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                color: white;
                padding: 40px;
                text-align: center;
                margin-bottom: 20px;
            }
            .app-header h1 {
                font-size: 2.5em;
                margin: 0 0 10px 0;
            }
            .app-header p {
                font-size: 1.1em;
                opacity: 0.9;
                margin: 0;
            }
            .stat-box {
                background: #e7f3ff;
                padding: 20px;
                border-radius: 8px;
                border-left: 4px solid #667eea;
                margin-bottom: 20px;
            }
            .stat-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
                gap: 15px;
            }
            .stat-item {
                background: white;
                padding: 15px;
                border-radius: 4px;
                text-align: center;
            }
            .stat-number {
                font-size: 1.8em;
                font-weight: bold;
                color: #667eea;
            }
            .stat-label {
                color: #666;
                font-size: 0.9em;
                margin-top: 5px;
            }
        """)
    ),
    
    ui.div(
        {"class": "app-header"},
        ui.h1("🧬 TA-like-ID Protein Analysis"),
        ui.p("Interactive filtering • C-terminal proximity analysis")
    ),
    
    ui.layout_sidebar(
        ui.sidebar(
            ui.h3("Filters"),
            
            ui.input_select(
                "operator",
                "C-terminal Distance",
                choices={
                    "<=": "≤ Less than or equal",
                    ">=": "≥ Greater than or equal",
                    "==": "== Exactly equal",
                    "<": "< Less than",
                    ">": "> Greater than"
                },
                selected="<="
            ),
            
            ui.input_numeric(
                "distance_value",
                "Threshold:",
                value=30,
                min=0,
                max=1000
            ),
            
            ui.input_select(
                "compare_column",
                "Compare by",
                choices={
                    "membrane_domain_count": "Membrane Domain Count",
                    "Prediction": "Prediction",
                },
                selected="membrane_domain_count"
            ),
            
            ui.input_action_button(
                "reset",
                "Reset Filters",
                class_="btn btn-secondary",
                style="width: 100%; margin-top: 10px;"
            ),
            
            ui.download_button(
                "download_csv",
                "📥 Download Filtered CSV",
                style="width: 100%; margin-top: 10px; background: #28a745; color: white;"
            )
        ),
        
        ui.output_ui("stats"),
        ui.output_plot("category_chart", height="400px"),
        ui.output_data_frame("data_table")
    )
)

# Server logic
def server(input: Inputs, output: Outputs, session: Session) -> None:
    
    @reactive.Calc
    def filtered_data() -> pd.DataFrame:
        """Apply filters to the data"""
        operator = input.operator()
        value = input.distance_value()
        
        if operator == "<=":
            mask = df['cterm_distance'] <= value
        elif operator == ">=":
            mask = df['cterm_distance'] >= value
        elif operator == "==":
            mask = df['cterm_distance'] == value
        elif operator == "<":
            mask = df['cterm_distance'] < value
        elif operator == ">":
            mask = df['cterm_distance'] > value
        else:
            mask = pd.Series([True] * len(df))
        
        return df[mask].copy()
    
    @reactive.Calc
    def category_counts() -> pd.Series:
        """Get category distribution"""
        fdf = filtered_data()
        compare_col = input.compare_column()
        
        if compare_col in fdf.columns:
            return fdf[compare_col].value_counts().head(10)
        return pd.Series(dtype=int)
    
    @output
    @render.ui
    def stats() -> Tag:
        """Render statistics summary"""
        total = len(df)
        filtered = len(filtered_data())
        reduction = ((1 - filtered/total) * 100) if total > 0 else 0
        
        return ui.div(
            {"class": "stat-box"},
            ui.h3("📊 Results Summary"),
            ui.div(
                {"class": "stat-grid"},
                ui.div(
                    {"class": "stat-item"},
                    ui.div({"class": "stat-number"}, f"{total:,}"),
                    ui.div({"class": "stat-label"}, "Total Proteins")
                ),
                ui.div(
                    {"class": "stat-item"},
                    ui.div({"class": "stat-number"}, f"{filtered:,}"),
                    ui.div({"class": "stat-label"}, "Filtered Results")
                ),
                ui.div(
                    {"class": "stat-item"},
                    ui.div({"class": "stat-number"}, f"{reduction:.1f}%"),
                    ui.div({"class": "stat-label"}, "Reduction")
                )
            )
        )
    
    @output
    @render.plot
    def category_chart() -> Figure:
        """Render category distribution chart"""
        counts = category_counts()
        
        if len(counts) == 0:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.text(0.5, 0.5, 'No data to display', 
                   ha='center', va='center', fontsize=14)
            ax.axis('off')
            return fig
        
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.bar(range(len(counts)), counts.values, color='#667eea', alpha=0.7, edgecolor='black')
        ax.set_xticks(range(len(counts)))
        ax.set_xticklabels(counts.index, rotation=45, ha='right')
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title(f'{input.compare_column()} Distribution', fontsize=14, fontweight='bold')
        ax.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        
        return fig
    
    @output
    @render.data_frame
    def data_table():
        """Render the filtered data table"""
        fdf = filtered_data()
        return render.DataGrid(fdf.head(100), width="100%", height="500px")
    
    @session.download(filename="filtered_proteins.csv")
    def download_csv():
        """Download filtered data as CSV"""
        fdf = filtered_data()
        csv_buffer = io.StringIO()
        fdf.to_csv(csv_buffer, index=False)
        csv_buffer.seek(0)
        return csv_buffer.getvalue()
    
    @reactive.Effect
    @reactive.event(input.reset)
    def _():
        """Reset filters to defaults"""
        ui.update_select("operator", selected="<=")
        ui.update_numeric("distance_value", value=30)
        ui.update_select("compare_column", selected="membrane_domain_count")

# Create the app
app = App(app_ui, server)
