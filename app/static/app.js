'use strict';

const byId = (id) => document.getElementById(id);
const ui = {
  form: byId('molecule-form'), input: byId('molecule-input'), error: byId('input-error'),
  generate: byId('generate-button'), generateLabel: byId('generate-label'),
  canvas: byId('structure-canvas'), svg: byId('svg-container'), placeholder: byId('canvas-placeholder'),
  hybridization: byId('hybridization-toggle'), groups: byId('functional-groups'),
  toggleGroups: byId('toggle-all-groups'), atomSelect: byId('atom-select'),
  status: byId('request-status'), warnings: byId('warnings'),
};
const examples = {
  aspirin: { input: 'CC(=O)Oc1ccccc1C(=O)O', format: 'smiles', name: 'Aspirin' },
  ethanol: { input: 'CH3CH2OH', format: 'condensed', name: 'Ethanol' },
  caffeine: { input: 'Cn1c(=O)c2c(ncn2C)n(C)c1=O', format: 'smiles', name: 'Caffeine' },
  isopropanol: { input: 'CH3CH(CH3)OH', format: 'condensed', name: 'Isopropanol' },
};
const state = {
  result: null, committed: null, groups: new Set(), selectedAtom: null,
  request: 0, controller: null, busy: false, svg: '',
};

function setStatus(message, error = false) {
  ui.status.textContent = message;
  ui.status.dataset.error = String(error);
}

function setError(message) {
  ui.error.textContent = message;
  ui.error.hidden = !message;
  ui.input.setAttribute('aria-invalid', message ? 'true' : 'false');
}

function showPlaceholder(message, loading = false) {
  ui.placeholder.replaceChildren();
  if (loading) {
    const orbit = document.createElement('div');
    orbit.className = 'loading-orbit';
    orbit.setAttribute('aria-hidden', 'true');
    ui.placeholder.append(orbit);
  }
  const text = document.createElement('p');
  text.textContent = message;
  ui.placeholder.append(text);
  ui.placeholder.hidden = false;
}

function setBusy(busy) {
  state.busy = busy;
  ui.canvas.setAttribute('aria-busy', String(busy));
  ui.generateLabel.textContent = busy ? 'Drawing your molecule…' : 'Generate structure';
  // A new input can supersede a request; drawing controls wait for its result.
  ui.hybridization.disabled = busy || !state.result;
  ui.atomSelect.disabled = busy || !state.result;
  ui.toggleGroups.disabled = busy || !state.result?.functional_groups.length;
  for (const button of ui.groups.querySelectorAll('button')) button.disabled = busy;
  byId('export-svg').disabled = busy || !state.result;
  byId('copy-smiles').disabled = busy || !state.result;
}

function currentFormat() {
  return ui.form.elements.format.value;
}

function setFormat(format) {
  const control = ui.form.querySelector(`input[name="format"][value="${format}"]`);
  if (control) control.checked = true;
}

function formulaFragment(formula) {
  const fragment = document.createDocumentFragment();
  for (const part of formula.match(/\d+|[^\d]+/g) || []) {
    if (/^\d+$/.test(part)) {
      const sub = document.createElement('sub');
      sub.textContent = part;
      fragment.append(sub);
    } else fragment.append(document.createTextNode(part));
  }
  return fragment;
}

function hybridizationLabel(value) {
  return ({ sp: 'sp', sp2: 'sp²', sp3: 'sp³', s: 's', other: 'Other' })[value] || value;
}

function setSelectedAtom(index) {
  const atom = state.result?.atoms.find((item) => item.index === index);
  state.selectedAtom = atom ? atom.index : null;
  ui.atomSelect.value = atom ? String(atom.index) : '';
  byId('selected-symbol').textContent = atom?.symbol || '·';
  byId('atom-hybridization').textContent = atom ? hybridizationLabel(atom.hybridization) : '—';
  byId('atom-hydrogens').textContent = atom ? String(atom.hydrogen_count) : '—';
  byId('atom-charge').textContent = atom ? (atom.formal_charge > 0 ? `+${atom.formal_charge}` : String(atom.formal_charge)) : '—';
  byId('atom-aromatic').textContent = atom ? (atom.aromatic ? 'Yes' : 'No') : '—';
  const explanations = {
    sp: 'Two directions of bonding; often associated with a linear environment.',
    sp2: 'Three directions of bonding; often associated with a planar environment.',
    sp3: 'Four directions of bonding; often associated with a tetrahedral environment.',
    s: 'An s-orbital assignment for this atom.',
    other: 'This atom is outside the common sp, sp², and sp³ descriptions.',
  };
  byId('hybridization-explanation').textContent = atom ? explanations[atom.hybridization] || 'Hybridization is a model of the local bonding environment.' : 'Choose an atom to see its bonding environment.';
  for (const selector of ui.svg.querySelectorAll('.atom-selector')) {
    selector.classList.toggle('selected-atom', Boolean(atom) && selector.classList.contains(`atom-${atom.index}`));
  }
}

function renderAtoms(result, reset) {
  ui.atomSelect.replaceChildren();
  const initial = document.createElement('option');
  initial.value = '';
  initial.textContent = 'Choose an atom in the drawing';
  ui.atomSelect.append(initial);
  for (const atom of result.atoms) {
    const option = document.createElement('option');
    option.value = String(atom.index);
    option.textContent = `Atom ${atom.index + 1} · ${atom.symbol} · ${hybridizationLabel(atom.hybridization)}`;
    ui.atomSelect.append(option);
  }
  ui.atomSelect.disabled = false;
  setSelectedAtom(reset ? null : state.selectedAtom);
}

function renderGroups(result) {
  const focusedGroup = ui.groups.contains(document.activeElement) ? document.activeElement.dataset.group : null;
  ui.groups.replaceChildren();
  byId('group-count').textContent = String(result.functional_groups.length);
  ui.toggleGroups.disabled = result.functional_groups.length === 0;
  const allHighlighted = result.functional_groups.length > 0 && result.functional_groups.every((group) => state.groups.has(group.id));
  ui.toggleGroups.textContent = allHighlighted ? 'Clear highlights' : 'Highlight all';
  if (!result.functional_groups.length) {
    const empty = document.createElement('p');
    empty.className = 'empty-text';
    empty.textContent = 'No supported functional groups found. The structure and atom details are still available.';
    ui.groups.append(empty);
  }
  for (const group of result.functional_groups) {
    const button = document.createElement('button');
    button.type = 'button';
    button.disabled = state.busy;
    button.className = 'group-button';
    button.dataset.group = group.id;
    button.setAttribute('aria-pressed', String(state.groups.has(group.id)));
    button.setAttribute('aria-label', `${group.name}, ${group.matches.length} ${group.matches.length === 1 ? 'match' : 'matches'}`);
    // Only accept hex colors from the chemistry service, never arbitrary CSS.
    const color = /^#[\da-f]{6}$/i.test(group.color) ? group.color : '#90a96d';
    button.style.setProperty('--group-color', color);
    button.style.setProperty('--group-tint', `${color}18`);
    const dot = document.createElement('span');
    dot.className = 'group-dot';
    dot.setAttribute('aria-hidden', 'true');
    const name = document.createElement('span');
    name.textContent = group.name;
    const count = document.createElement('span');
    count.className = 'group-matches';
    count.textContent = String(group.matches.length);
    const check = document.createElement('span');
    check.className = 'group-check';
    check.setAttribute('aria-hidden', 'true');
    check.textContent = state.groups.has(group.id) ? '' : '+';
    check.classList.toggle('is-checked', state.groups.has(group.id));
    button.append(dot, name, count, check);
    ui.groups.append(button);
    if (focusedGroup === group.id) button.focus({ preventScroll: true });
  }
}

function renderResult(result, request, reset) {
  state.result = result;
  state.committed = request;
  const validGroups = new Set(result.functional_groups.map((group) => group.id));
  state.groups = new Set([...state.groups].filter((id) => validGroups.has(id)));
  state.svg = result.svg;
  // SVG comes only from our same-origin RDKit drawing endpoint, never from input markup.
  ui.svg.innerHTML = result.svg;
  const svg = ui.svg.querySelector('svg');
  if (svg) {
    svg.removeAttribute('width');
    svg.removeAttribute('height');
    svg.setAttribute('role', 'img');
    svg.setAttribute('aria-label', `Skeletal structure of ${result.molecular_formula}. Choose an atom using the atom inspector below.`);
  }
  ui.placeholder.hidden = true;
  byId('structure-title').textContent = request.name || 'Your molecule';
  byId('molecule-subtitle').textContent = `${result.input_format === 'condensed' ? 'Condensed formula' : 'SMILES'} → skeletal structure`;
  byId('metric-formula').replaceChildren(formulaFragment(result.molecular_formula));
  byId('metric-mass').textContent = Number(result.molecular_weight).toFixed(2);
  byId('metric-atoms').textContent = String(result.atom_count);
  byId('metric-bonds').textContent = String(result.bond_count);
  byId('canonical-smiles').textContent = result.canonical_smiles;
  byId('canonical-smiles').title = result.canonical_smiles;
  ui.warnings.replaceChildren();
  for (const warning of result.warnings || []) {
    const item = document.createElement('p');
    item.textContent = warning;
    ui.warnings.append(item);
  }
  ui.warnings.hidden = !result.warnings?.length;
  renderGroups(result);
  renderAtoms(result, reset);
}

async function analyze(request, reset = true) {
  const input = request.input.trim();
  if (!input) {
    setError('Enter a SMILES string or a condensed structural formula to begin.');
    ui.input.focus();
    return;
  }
  if (reset) {
    state.groups = new Set();
    setError('');
  }
  state.controller?.abort();
  state.controller = new AbortController();
  const currentRequest = ++state.request;
  const chosenGroups = [...state.groups];
  const showHybridization = ui.hybridization.checked;
  const focusedElement = document.activeElement;
  const focusedGroup = focusedElement?.dataset.group;
  const focusedId = focusedElement?.id;
  setBusy(true);
  showPlaceholder('Drawing your molecule…', true);
  setStatus(reset ? 'Building the structure and checking its chemistry…' : 'Updating the structure…');
  try {
    const response = await fetch('/api/analyze', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ input, format: request.format, show_hybridization: showHybridization, highlight_groups: chosenGroups }),
      signal: state.controller.signal,
    });
    let result;
    try { result = await response.json(); }
    catch { throw new Error('The drawing service returned an unexpected response. Please try again.'); }
    if (!response.ok) {
      throw new Error(typeof result.detail === 'string' ? result.detail : 'We couldn’t read that structure. Check the input and its selected format.');
    }
    if (currentRequest !== state.request) return;
    renderResult(result, { ...request, input }, reset);
    const count = result.functional_groups.length;
    setStatus(`Structure ready · ${result.atom_count} heavy atoms · ${count} functional group ${count === 1 ? 'type' : 'types'} found.`);
  } catch (error) {
    if (error.name === 'AbortError' || currentRequest !== state.request) return;
    const message = error instanceof TypeError ? 'Couldn’t reach the drawing service. Check your connection and try again.' : error.message;
    if (reset) setError(message);
    setStatus(`${message}${state.result ? ' The previous molecule is still shown.' : ''}`, true);
    if (state.result) {
      ui.placeholder.hidden = true;
      // Put controls back in sync with the last successful drawing, including
      // when a new input failed after the previous molecule had highlights.
      ui.hybridization.checked = state.renderedHybridization;
      state.groups = new Set(state.renderedGroups || []);
      renderGroups(state.result);
    } else showPlaceholder('Your structure will appear here. Check the input and try again.');
  } finally {
    if (currentRequest === state.request) {
      setBusy(false);
      if (document.activeElement === document.body) {
        const focusTarget = focusedGroup
          ? [...ui.groups.querySelectorAll('button')].find((button) => button.dataset.group === focusedGroup)
          : focusedId ? byId(focusedId) : null;
        focusTarget?.focus({ preventScroll: true });
      }
      if (!ui.status.dataset.error || ui.status.dataset.error === 'false') {
        state.renderedHybridization = showHybridization;
        state.renderedGroups = chosenGroups;
      }
    }
  }
}

ui.form.addEventListener('submit', (event) => {
  event.preventDefault();
  analyze({ input: ui.input.value, format: currentFormat(), name: null });
});
ui.input.addEventListener('keydown', (event) => {
  if (event.key === 'Enter' && !event.shiftKey && !event.isComposing) {
    event.preventDefault();
    ui.form.requestSubmit();
  }
});
ui.input.addEventListener('input', () => {
  setError('');
  if (state.committed && ui.input.value.trim() !== state.committed.input) {
    setStatus('Input changed. Generate a structure to update the drawing.');
  }
});
ui.form.querySelectorAll('input[name="format"]').forEach((radio) => radio.addEventListener('change', () => {
  setError('');
  if (state.committed) setStatus('Format changed. Generate a structure to update the drawing.');
}));
byId('clear-input').addEventListener('click', () => {
  ui.input.value = '';
  setError('');
  setStatus('Enter a new molecule. The current drawing stays available until you generate another.');
  ui.input.focus();
});
document.querySelectorAll('[data-example]').forEach((button) => button.addEventListener('click', () => {
  const example = examples[button.dataset.example];
  ui.input.value = example.input;
  setFormat(example.format);
  analyze(example);
}));
ui.hybridization.addEventListener('change', () => {
  if (state.committed) analyze(state.committed, false);
});
ui.groups.addEventListener('click', (event) => {
  const button = event.target.closest('button[data-group]');
  if (!button || !state.committed || state.busy) return;
  const group = button.dataset.group;
  state.groups.has(group) ? state.groups.delete(group) : state.groups.add(group);
  renderGroups(state.result);
  analyze(state.committed, false);
});
ui.toggleGroups.addEventListener('click', () => {
  if (!state.committed || !state.result || state.busy) return;
  const all = state.result.functional_groups.map((group) => group.id);
  state.groups = all.every((id) => state.groups.has(id)) ? new Set() : new Set(all);
  renderGroups(state.result);
  analyze(state.committed, false);
});
ui.atomSelect.addEventListener('change', () => setSelectedAtom(ui.atomSelect.value === '' ? null : Number(ui.atomSelect.value)));
ui.svg.addEventListener('click', (event) => {
  let target = event.target;
  while (target && target !== ui.svg) {
    const match = target.getAttribute?.('class')?.match(/(?:^|\s)atom-(\d+)(?:\s|$)/);
    // Bond paths contain two atom classes; only atom labels and selectors select an atom.
    if (match && !/(?:^|\s)bond-\d+/.test(target.getAttribute('class') || '')) {
      setSelectedAtom(Number(match[1]));
      setStatus(`Selected atom ${Number(match[1]) + 1}. Details are in the atom inspector.`);
      break;
    }
    target = target.parentElement;
  }
});
byId('export-svg').addEventListener('click', () => {
  if (!state.svg) return;
  const url = URL.createObjectURL(new Blob([state.svg], { type: 'image/svg+xml;charset=utf-8' }));
  const link = document.createElement('a');
  link.href = url;
  const basename = (state.committed?.name || state.result.molecular_formula || 'molecule').toLowerCase().replace(/[^a-z0-9_-]+/g, '-');
  link.download = `${basename}-structure.svg`;
  document.body.append(link);
  link.click();
  link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
  setStatus('Saved the current structure as an SVG image.');
});
byId('copy-smiles').addEventListener('click', async () => {
  if (!state.result) return;
  try {
    await navigator.clipboard.writeText(state.result.canonical_smiles);
    setStatus('Canonical SMILES copied to your clipboard.');
  } catch {
    setStatus('Clipboard access is unavailable. Select the canonical SMILES text to copy it.');
  }
});

analyze(examples.aspirin);
