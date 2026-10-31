"""HTTP interface and static frontend for carbo cat."""

from pathlib import Path
from typing import Literal

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, ConfigDict, Field

from app.chemistry import analyze_molecule

STATIC = Path(__file__).resolve().parent / "static"

app = FastAPI(
    title="carbo cat",
    description="Draw molecules, inspect atom hybridization, and find functional groups.",
    version="0.1.0",
)


class AnalyzeRequest(BaseModel):
    model_config = ConfigDict(extra="forbid", str_strip_whitespace=True)

    input: str = Field(min_length=1, max_length=500)
    format: Literal["auto", "smiles", "condensed"] = "auto"
    show_hybridization: bool = False
    highlight_groups: list[str] = Field(default_factory=list, max_length=30)


@app.get("/api/health")
def health():
    return {"status": "ok", "version": "0.1.0"}


@app.post("/api/analyze")
def analyze(request: AnalyzeRequest):
    """Analyze explicit connectivity; return an SVG plus chemical metadata."""
    try:
        return analyze_molecule(
            request.input,
            input_format=request.format,
            show_hybridization=request.show_hybridization,
            highlight_groups=request.highlight_groups,
        )
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.get("/", include_in_schema=False)
def index():
    return FileResponse(STATIC / "index.html")


app.mount("/static", StaticFiles(directory=STATIC), name="static")
