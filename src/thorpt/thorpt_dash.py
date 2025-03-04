import dash
from dash import dcc, html, Input, Output, State
import dash_bootstrap_components as dbc
import pandas as pd
import subprocess
import os
import webbrowser
from threading import Timer

def open_browser():
    webbrowser.open_new("http://127.0.0.1:8050")

# Initialize Dash app with Bootstrap theme
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Layout
app.layout = dbc.Container([
    html.H1("Python Script Runner", className="text-center mt-4 mb-2"),
    
    # File upload
    dcc.Upload(
        id='upload-file',
        children=dbc.Button("Upload TXT File", color="primary", className="mb-3"),
        multiple=False
    ),
    html.Div(id="file-status", className="mb-3"),

    # Parameter input
    dbc.Input(id="param-input", type="text", placeholder="Enter parameters...", className="mb-3"),

    # Run button
    dbc.Button("Run Script", id="run-btn", color="success", className="mb-3"),
    
    # Output display
    dbc.Alert(id="output-display", color="info", className="mt-3", is_open=False)
], fluid=True)


# Callback for file upload
@app.callback(
    Output("file-status", "children"),
    Input("upload-file", "contents"),
    State("upload-file", "filename")
)
def save_uploaded_file(contents, filename):
    if contents:
        filepath = os.path.join(os.getcwd(), filename)
        with open(filepath, "wb") as f:
            f.write(contents.encode('utf-8'))  # Save file
        return f"File uploaded: {filename}"
    return ""


# Callback for running the script
@app.callback(
    Output("output-display", "children"),
    Output("output-display", "is_open"),
    Input("run-btn", "n_clicks"),
    State("param-input", "value"),
    prevent_initial_call=True
)
def run_script(n_clicks, params):
    try:
        # Modify this command to run your specific script
        result = subprocess.run(["python", "your_script.py", "--from-dash", params], capture_output=True, text=True)
        return f"Script Output: {result.stdout}", True
    except Exception as e:
        return f"Error: {str(e)}", True


if __name__ == '__main__':
    Timer(1, open_browser).start()  # Open browser after 1 second
    app.run_server(debug=False, host='127.0.0.1', port=8050)
