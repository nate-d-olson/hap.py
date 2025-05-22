#!/usr/bin/env python3
"""
Generate a visual dashboard for the hap.py modernization progress.

This script creates an HTML dashboard showing the current status of 
the modernization effort, including code quality metrics, C++ dependency
reduction, and test coverage.
"""

import argparse
import os
from datetime import datetime
from pathlib import Path

# Try to import track_progress module
try:
    from track_progress import ModernizationTracker
except ImportError:
    # Fall back to assuming it's in the same directory
    import sys
    sys.path.append(str(Path(__file__).parent))
    from track_progress import ModernizationTracker


class DashboardGenerator:
    """Generate a visual dashboard for modernization progress."""
    
    def __init__(self, project_root: Path):
        """
        Initialize the dashboard generator.
        
        Args:
            project_root: Root directory of the project
        """
        self.project_root = project_root
        self.tracker = ModernizationTracker()
        self.dashboard_dir = project_root / "dashboard"
        self.ensure_dashboard_dir()
    
    def ensure_dashboard_dir(self):
        """Ensure the dashboard directory exists."""
        os.makedirs(self.dashboard_dir, exist_ok=True)
    
    def generate_dashboard(self, output_path: str = None) -> str:
        """
        Generate the dashboard HTML.
        
        Args:
            output_path: Path to save the HTML file
            
        Returns:
            Path to the generated HTML file
        """
        # Update metrics
        self.tracker.update_metrics()
        
        # Generate HTML
        html = self.generate_html()
        
        # Save HTML
        if output_path is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_path = self.dashboard_dir / f"dashboard_{timestamp}.html"
        else:
            output_path = Path(output_path)
        
        with open(output_path, "w") as f:
            f.write(html)
        
        print(f"✅ Dashboard generated: {output_path}")
        return str(output_path)
    
    def generate_html(self) -> str:
        """
        Generate the dashboard HTML.
        
        Returns:
            HTML content
        """
        # Get metrics
        quality = self.tracker.check_code_quality()
        cpp_deps = self.tracker.check_cpp_dependencies()
        packages = self.tracker.check_dependency_packages()
        
        # Calculate overall progress
        total_cpp_components = len(cpp_deps)
        migrated_components = sum(1 for comp in cpp_deps.values() if comp["has_python_replacement"])
        migration_percentage = (migrated_components / total_cpp_components * 100) if total_cpp_components > 0 else 0
        
        # Get Python 2 artifacts
        artifacts = quality["python2_artifacts"]
        total_artifacts = sum(artifacts.values())
        
        # Generate HTML
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>hap.py Modernization Dashboard</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
            margin: 0;
            padding: 20px;
            color: #333;
        }}
        h1, h2, h3 {{
            color: #2c3e50;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
        }}
        .dashboard {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
        }}
        .card {{
            background: #fff;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            padding: 20px;
            margin-bottom: 20px;
        }}
        .progress-bar {{
            background-color: #f1f1f1;
            border-radius: 4px;
            height: 24px;
            margin: 10px 0;
            position: relative;
        }}
        .progress-bar .fill {{
            background-color: #4caf50;
            height: 100%;
            border-radius: 4px;
            transition: width 0.3s ease;
        }}
        .progress-bar .label {{
            position: absolute;
            top: 0;
            left: 0;
            right: 0;
            line-height: 24px;
            text-align: center;
            font-weight: bold;
            color: #333;
        }}
        .status {{
            padding: 4px 8px;
            border-radius: 4px;
            display: inline-block;
            margin-left: 5px;
        }}
        .success {{
            background-color: #dff0d8;
            color: #3c763d;
        }}
        .warning {{
            background-color: #fcf8e3;
            color: #8a6d3b;
        }}
        .danger {{
            background-color: #f2dede;
            color: #a94442;
        }}
        .info {{
            background-color: #d9edf7;
            color: #31708f;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
        }}
        table th, table td {{
            padding: 8px 12px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }}
        table th {{
            background-color: #f5f5f5;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>hap.py Modernization Dashboard</h1>
        <p>Updated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
        
        <div class="card">
            <h2>Overall Progress</h2>
            <div class="progress-bar">
                <div class="fill" style="width: {migration_percentage}%;"></div>
                <div class="label">{migration_percentage:.1f}% ({migrated_components}/{total_cpp_components} components)</div>
            </div>
        </div>
        
        <div class="dashboard">
            <!-- Left Column -->
            <div>
                <div class="card">
                    <h2>Code Quality Status</h2>
                    <p>
                        <strong>Pre-commit:</strong> 
                        {('<span class="status success">✓ Configured</span>' if quality['pre_commit_installed'] else '<span class="status danger">✗ Not configured</span>')}
                    </p>
                    <p>
                        <strong>Type Coverage:</strong> 
                        <div class="progress-bar">
                            <div class="fill" style="width: {quality['type_coverage']}%;"></div>
                            <div class="label">{quality['type_coverage']}%</div>
                        </div>
                    </p>
                    <p>
                        <strong>Python 2 Artifacts:</strong> 
                        {('<span class="status success">✓ None found</span>' if total_artifacts == 0 else f'<span class="status danger">✗ {total_artifacts} issues found</span>')}
                    </p>
                    {self._generate_artifacts_table(artifacts) if total_artifacts > 0 else ''}
                </div>
                
                <div class="card">
                    <h2>Modern Package Adoption</h2>
                    <table>
                        <tr>
                            <th>Package</th>
                            <th>Status</th>
                        </tr>
                        {self._generate_package_rows(packages)}
                    </table>
                </div>
            </div>
            
            <!-- Right Column -->
            <div>
                <div class="card">
                    <h2>C++ Component Migration</h2>
                    <table>
                        <tr>
                            <th>Component</th>
                            <th>Status</th>
                            <th>Files</th>
                        </tr>
                        {self._generate_component_rows(cpp_deps)}
                    </table>
                </div>
                
                <div class="card">
                    <h2>Testing Status</h2>
                    <p><strong>Pytest Files:</strong> {quality['test_files']['pytest_files']}</p>
                    <p><strong>Shell Test Files:</strong> {quality['test_files']['shell_test_files']}</p>
                    <p><strong>Total Test Files:</strong> {quality['test_files']['total_tests']}</p>
                </div>
                
                <div class="card">
                    <h2>Next Priority Actions</h2>
                    <ol>
                        {self._generate_priority_actions(quality, cpp_deps, packages, total_artifacts)}
                    </ol>
                </div>
            </div>
        </div>
    </div>
</body>
</html>
"""
        return html
    
    def _generate_artifacts_table(self, artifacts):
        """Generate HTML table for Python 2 artifacts."""
        rows = ""
        for artifact, count in artifacts.items():
            if count > 0:
                status_class = "danger" if count > 10 else "warning"
                rows += f"""
                <tr>
                    <td>{artifact.replace("_", " ").title()}</td>
                    <td><span class="status {status_class}">{count}</span></td>
                </tr>
                """
        
        if rows:
            return f"""
            <table>
                <tr>
                    <th>Artifact Type</th>
                    <th>Count</th>
                </tr>
                {rows}
            </table>
            """
        else:
            return ""
    
    def _generate_package_rows(self, packages):
        """Generate HTML table rows for package adoption."""
        rows = ""
        for package, adopted in packages.items():
            status = '<span class="status success">✓ Adopted</span>' if adopted else '<span class="status danger">✗ Missing</span>'
            rows += f"""
            <tr>
                <td>{package}</td>
                <td>{status}</td>
            </tr>
            """
        return rows
    
    def _generate_component_rows(self, cpp_deps):
        """Generate HTML table rows for C++ components."""
        rows = ""
        for component, details in cpp_deps.items():
            if details["has_python_replacement"]:
                status = '<span class="status success">✓ Python replacement</span>'
                status_class = "success"
            elif details["total_files"] == 0:
                status = '<span class="status info">ℹ Not found</span>'
                status_class = "info"
            else:
                status = '<span class="status warning">⟳ Needs migration</span>'
                status_class = "warning"
            
            rows += f"""
            <tr>
                <td>{component}</td>
                <td>{status}</td>
                <td>{details["total_files"]}</td>
            </tr>
            """
        return rows
    
    def _generate_priority_actions(self, quality, cpp_deps, packages, total_artifacts):
        """Generate HTML list items for priority actions."""
        actions = []
        
        if total_artifacts > 0:
            actions.append('<li><strong>Fix Python 2 artifacts</strong> - Run <code>pre-commit run pyupgrade --all-files</code></li>')
        
        if not quality['pre_commit_installed']:
            actions.append('<li><strong>Set up pre-commit</strong> - Run <code>pip install pre-commit && pre-commit install</code></li>')
        
        # Find components that should be prioritized
        high_priority_components = []
        for component, details in cpp_deps.items():
            if not details["has_python_replacement"] and details["total_files"] > 0:
                # Prioritize by component importance
                if component in ["blocksplit", "quantify", "vcfcheck", "preprocess"]:
                    high_priority_components.append(component)
        
        if high_priority_components:
            components_list = ", ".join(high_priority_components)
            actions.append(f'<li><strong>Migrate high-priority components</strong> - Focus on: {components_list}</li>')
        
        # Check for important packages
        missing_packages = []
        for package in ["pysam", "biopython", "pytest"]:
            if not packages.get(package, False):
                missing_packages.append(package)
        
        if missing_packages:
            packages_list = ", ".join(missing_packages)
            actions.append(f'<li><strong>Add missing dependencies</strong> - Install: {packages_list}</li>')
        
        return "\n".join(actions) or "<li>No immediate actions - keep up the good work!</li>"


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Generate a visual dashboard for hap.py modernization progress"
    )
    parser.add_argument(
        "-o", "--output", 
        help="Path to save the HTML file (default: dashboard/dashboard_TIMESTAMP.html)"
    )
    parser.add_argument(
        "--open", action="store_true",
        help="Open the dashboard in a web browser after generation"
    )
    
    args = parser.parse_args()
    
    # Get project root directory (parent of the script directory)
    project_root = Path(__file__).parent.parent
    
    # Generate dashboard
    generator = DashboardGenerator(project_root)
    dashboard_path = generator.generate_dashboard(args.output)
    
    # Open in browser if requested
    if args.open:
        try:
            import webbrowser
            webbrowser.open(f"file://{os.path.abspath(dashboard_path)}")
        except Exception as e:
            print(f"Could not open browser: {e}")


if __name__ == "__main__":
    main()
