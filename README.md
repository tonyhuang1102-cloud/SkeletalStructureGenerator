# carbo cat

A cat-themed chemistry website for turning **SMILES** or supported **condensed structural formulas** into skeletal structures. Explore atom hybridization, highlight functional groups, and download a structure as SVG.

Built with Python, RDKit, and FastAPI, with a plain HTML/CSS/JavaScript frontend. No API keys, database, Node build, or external chemistry service are required. Molecule analysis happens on the server; the application does not save submitted molecules.

## Open it as a public website

[Deploy carbo cat on Render](https://render.com/deploy?repo=https://github.com/tonyhuang1102-cloud/SkeletalStructureGenerator)

The link opens Render's deployment flow for this repository. Sign in with GitHub, authorize access to the repository if prompted, review the free service configuration, and deploy. The website address is assigned by Render after the service is created.

The source files must first be saved to your GitHub repository. A development server inside a Codex workspace does not provide a public website address.

1. Upload the project contents to your GitHub repository. Keep `render.yaml`, `requirements.txt`, and the `app` folder at the repository's top level, rather than inside an extra folder.
2. Sign in to [Render](https://dashboard.render.com) using GitHub.
3. Choose **New → Blueprint**, connect this repository, and let Render read `render.yaml`.
4. Review the proposed service and plan, then deploy. No API keys or secret environment variables are needed.
5. When the service reports that it is live, open the website address Render displays. That is your shareable carbo cat link.

The Blueprint requests a free service. Check Render's current availability and terms before confirming deployment; free instances can sleep when idle, so the first visit may take longer. Deployment has not been performed or verified merely by adding this configuration.

## Run in the prepared cloud environment

```bash
cd /workspace/SkeletalStructureGenerator
source /workspace/.venvs/skeletal/bin/activate
python -m uvicorn app.main:app --host 0.0.0.0 --port 8000
```

The application listens on port **8000**. On your own computer, open `http://localhost:8000`. A cloud workspace requires its port preview or forwarding feature; a public website needs a separate deployment.

## Set up on another computer

Use Python **3.12** and a platform supported by RDKit's Python wheels.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-dev.txt
python -m uvicorn app.main:app --reload
```

On Windows, activate with `.venv\Scripts\activate` instead. For running without development tests, install `requirements.txt`.

## Use it

1. Choose Auto, SMILES, or Condensed, enter a structure, and select **Generate**. Example buttons provide molecules to explore.
2. Turn on hybridization labels to see RDKit's atom assignments. Select an atom in the drawing or atom list for its details.
3. Select functional groups to highlight matching atoms and bonds.
4. Download an SVG to keep the displayed structure.

Try:

| Molecule | Format | Input |
| --- | --- | --- |
| Ethanol | Condensed | `CH3CH2OH` |
| Acetic acid | Condensed | `CH3COOH` |
| Acetone | Condensed | `CH3COCH3` |
| Isopropanol | Condensed | `CH3CH(CH3)OH` |
| Aspirin | SMILES | `CC(=O)Oc1ccccc1C(=O)O` |
| Benzene | SMILES | `c1ccccc1` |
| Acetonitrile | SMILES | `CC#N` |

### Chemistry scope

- SMILES supports explicit molecular connectivity, including aromatic rings, bond orders, charges, and stereochemical input handled by RDKit. The output is a 2D depiction, not a 3D geometry prediction.
- Condensed input is a bounded subset of structural notation with explicit atom order, simple branches, repeated groups, and common carbonyl notation. It is not a universal chemical-name or formula resolver.
- A molecular formula such as `C2H6O` does not specify whether the molecule is ethanol or dimethyl ether. Enter a structural formula or SMILES instead.
- Hybridization is RDKit's model assignment. It is useful for learning and inspection, but it does not represent a measured orbital structure. Some unusual atoms do not have an sp/sp²/sp³ assignment.
- Functional groups are identified by a curated set of SMARTS patterns. Unlisted groups are not identified; matching covers the supported patterns, not every chemical motif.
- Input and molecule-size limits keep interactive analysis small. Very large structures and excessive group expansion are rejected with a validation error.

## Tests

From the repository root, with the Python environment activated:

```bash
python -m pytest -q
```

Tests cover condensed parsing, invalid/ambiguous structures, chemistry assignments, functional-group exclusions, rendering, and the HTTP API. Interactive browser checks can additionally exercise input, toggles, atom selection, SVG download, and the mobile layout.

## API

`GET /api/health` provides a readiness check. `POST /api/analyze` accepts:

```json
{
  "input": "CC(=O)O",
  "format": "smiles",
  "show_hybridization": true,
  "highlight_groups": []
}
```

The response includes canonical SMILES, molecular formula and mass, atom and bond counts, per-atom properties, functional-group matches, warnings, and an SVG. Use group IDs from the response in `highlight_groups`. Format is `auto`, `smiles`, or `condensed`. Invalid inputs return HTTP 422. Interactive API documentation is at `/docs` while the app is running.

## Project files

- `app/condensed.py`: condensed structural formula parser.
- `app/chemistry.py`: RDKit analysis, group matching, and SVG drawing.
- `app/main.py`: API and static-file serving.
- `app/static/`: frontend.
- `tests/`: chemistry and API regression tests.

This is a first development version. Public hosting, authentication, and production traffic controls are separate deployment work. Publishing the Codex development environment does not publish this website.
