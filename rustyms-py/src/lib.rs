//! Python bindings to the rustyms library.

use std::fmt::Debug;

use ordered_float::OrderedFloat;
use pyo3::prelude::*;
use pyo3::{exceptions::PyValueError, types::PyType};

use rustyms::{Chemical, MultiChemical};

/// Mass mode enum.
#[pyclass]
enum MassMode {
    Monoisotopic,
    Average,
}

/// Element.
///
/// A chemical element, with its isotopes and their properties.
///
/// Parameters
/// ----------
/// symbol : str
///
#[pyclass]
#[derive(Debug, Clone)]
pub struct Element(rustyms::Element);

#[pymethods]
impl Element {
    #[new]
    fn new(symbol: &str) -> PyResult<Self> {
        rustyms::Element::try_from(symbol)
            .map(Element)
            .map_err(|_| PyValueError::new_err("Invalid element symbol."))
    }

    fn __repr__(&self) -> String {
        format!("Element('{}')", self.0)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    /// Get all available isotopes (N, mass, abundance).
    ///
    /// Returns
    /// -------
    /// list[tuple[int, float, float]]
    ///
    fn isotopes(&self) -> Vec<(u16, f64, f64)> {
        self.0
            .isotopes()
            .iter()
            .map(|i| (i.0, i.1.value, i.2))
            .collect()
    }

    /// The mass of the specified isotope of this element (if that isotope exists).
    ///
    /// Parameters
    /// ----------
    /// isotope : int | None
    ///    The isotope number (default: None).
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    fn mass(&self, isotope: Option<u16>) -> Option<f64> {
        self.0.mass(isotope).map(|mass| mass.value)
    }

    /// The average weight of the specified isotope of this element (if that isotope exists).
    ///
    /// Parameters
    /// ----------
    /// isotope : int | None
    ///     The isotope number (default: None).
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn average_weight(&self, isotope: Option<u16>) -> Option<f64> {
        self.0.average_weight(isotope).map(|mass| mass.value)
    }
}

impl std::fmt::Display for Element {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Molecular formula.
///
/// A molecular formula: a selection of elements of specified isotopes together forming a structure.
///
#[pyclass]
#[derive(Debug, Clone)]
pub struct MolecularFormula(rustyms::MolecularFormula);

#[pymethods]
impl MolecularFormula {
    /// Create a new molecular formula.
    ///
    /// Parameters
    /// ----------
    /// elements : list[tuple[Element, int | None, int]]
    ///     List of tuples of elements, isotope numbers and counts.
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    /// Raises
    /// ------
    /// ValueError
    ///    If an isotope is invalid.
    ///
    #[new]
    fn new(elements: Vec<(Element, Option<u16>, i16)>) -> PyResult<Self> {
        let elements = elements
            .iter()
            .map(|(e, i, n)| (e.0, *i, *n))
            .collect::<Vec<_>>();
        let formula = rustyms::MolecularFormula::new(elements.as_slice())
            .ok_or_else(|| PyValueError::new_err("Invalid isotopes"))?;
        Ok(MolecularFormula(formula))
    }

    fn __repr__(&self) -> String {
        format!("MolecularFormula({})", self.0)
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    /// Create a new molecular formula from a ProForma formula notation string.
    ///
    /// Parameters
    /// ----------
    /// proforma : str
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    #[classmethod]
    fn from_pro_forma(_cls: &PyType, proforma: &str) -> PyResult<Self> {
        rustyms::MolecularFormula::from_pro_forma(proforma)
            .map(MolecularFormula)
            .map_err(|e| PyValueError::new_err(format!("Invalid ProForma string: {}", e)))
    }

    /// Create a new molecular formula from a PSI-MOD formula notation string.
    ///
    /// Parameters
    /// ----------
    /// psi_mod : str
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    #[classmethod]
    fn from_psi_mod(_cls: &PyType, psi_mod: &str) -> PyResult<Self> {
        rustyms::MolecularFormula::from_psi_mod(psi_mod)
            .map(MolecularFormula)
            .map_err(|e| PyValueError::new_err(format!("Invalid PSI-MOD string: {}", e)))
    }

    /// Add the given element to this formula (while keeping it ordered and simplified)
    ///
    /// Parameters
    /// ----------
    /// element : Element
    ///     The element to add.
    /// isotope : int | None
    ///     The isotope number of the element to add (default: None).
    /// n : int
    ///     The number of atoms of this element to add (default: 1).
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the element or isotope is invalid.
    ///
    #[pyo3(signature = (element, isotope=None, n=1))]
    fn add(&mut self, element: &Element, isotope: Option<u16>, n: i16) -> PyResult<Option<()>> {
        if self.0.add((element.0, isotope, n)) {
            Ok(None)
        } else {
            Err(PyValueError::new_err("Invalid element or isotope"))
        }
    }

    /// Get the elements making this formula.
    ///
    /// Returns
    /// -------
    /// list[tuple[Element, int | None, int]]
    ///
    fn elements(&self) -> Vec<(Element, Option<u16>, i16)> {
        self.0
            .elements()
            .iter()
            .map(|(e, i, n)| (Element(*e), *i, *n))
            .collect()
    }

    /// Create a new molecular formula with the given global isotope modifications.
    fn with_global_isotope_modifications(
        &self,
        substitutions: Vec<(Element, Option<u16>)>,
    ) -> PyResult<Self> {
        match self.0.with_global_isotope_modifications(
            substitutions
                .iter()
                .map(|(e, i)| (e.0, *i))
                .collect::<Vec<_>>()
                .as_slice(),
        ) {
            Some(formula) => Ok(MolecularFormula(formula)),
            None => Err(PyValueError::new_err(
                "Invalid global isotope modifications",
            )),
        }
    }

    /// Get the number of electrons (the only charged species, any ionic species is saved as that element +/- the correct number of electrons). The inverse of that number is given as the charge.
    ///
    /// Returns
    /// -------
    /// int
    ///
    fn charge(&self) -> i16 {
        self.0.charge()
    }

    /// The mass of the molecular formula of this element, if all element species (isotopes) exists
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn monoisotopic_mass(&self) -> f64 {
        self.0.monoisotopic_mass().value
    }

    /// The average weight of the molecular formula of this element, if all element species (isotopes) exists.
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn average_weight(&self) -> f64 {
        self.0.average_weight().value
    }

    /// Get the mass in the given mode.
    ///
    /// Parameters
    /// ----------
    /// mode : MassMode
    ///    The mode to get the mass in.
    ///
    /// Returns
    /// -------
    /// float
    ///
    /// Raises
    /// ------
    /// ValueError
    ///   If the mode is not one of the valid modes.
    ///
    #[pyo3(signature = (mode=&MassMode::Monoisotopic))]
    fn mass(&self, mode: &MassMode) -> PyResult<f64> {
        match mode {
            MassMode::Monoisotopic => Ok(self.monoisotopic_mass()),
            MassMode::Average => Ok(self.average_weight()),
        }
    }

    /// Create a Hill notation from this collections of elements merged with the pro forma notation for specific isotopes.
    ///
    /// Returns
    /// -------
    /// str
    fn hill_notation(&self) -> String {
        self.0.hill_notation()
    }

    /// Create a Hill notation from this collections of elements merged with the pro forma notation for specific isotopes. Using fancy unicode characters for subscript and superscript numbers.
    ///
    /// Returns
    /// -------
    /// str
    fn hill_notation_fancy(&self) -> String {
        self.0.hill_notation_fancy()
    }

    /// Create a Hill notation from this collections of elements encoded in HTML.
    ///
    /// Returns
    /// -------
    /// str
    fn hill_notation_html(&self) -> String {
        self.0.hill_notation_html()
    }
}

/// A selection of ions that together define the charge of a peptide.
#[pyclass]
pub struct MolecularCharge(rustyms::MolecularCharge);

#[pymethods]
impl MolecularCharge {
    /// Create a charge state with the given ions.
    ///
    /// Parameters
    /// ----------
    /// charge_carriers : list[tuple[int, MolecularFormula]]
    ///   The charge carriers.
    ///
    /// Returns
    /// -------
    /// MolecularCharge
    ///
    #[new]
    fn new(charge_carriers: Vec<(i32, MolecularFormula)>) -> Self {
        MolecularCharge(rustyms::MolecularCharge {
            charge_carriers: charge_carriers
                .iter()
                .map(|(n, mol)| (*n as isize, mol.0.clone()))
                .collect(),
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "MolecularCharge(charge_carriers={})",
            self.0
                .charge_carriers
                .iter()
                .map(|(n, mol)| format!("({}, {})", n, mol))
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    // Create a default charge state with only protons.
    ///
    /// Parameters
    /// ----------
    /// charge : int
    ///    The charge.
    ///
    /// Returns
    /// -------
    /// MolecularCharge
    ///
    #[classmethod]
    fn proton(_cls: &PyType, charge: i32) -> Self {
        MolecularCharge(rustyms::MolecularCharge::proton(charge as isize))
    }

    /// List of counts and molecular formulas for the charge carriers.
    ///
    /// Returns
    /// -------
    /// list[tuple[int, MolecularFormula]]
    ///     The charge carriers.
    #[getter]
    fn charge_carriers(&self) -> Vec<(i32, MolecularFormula)> {
        self.0
            .charge_carriers
            .iter()
            .map(|(n, mol)| (*n as i32, MolecularFormula(mol.clone())))
            .collect()
    }
}

/// Amino acid.
///
/// Parameters
/// ----------
/// name : str
///    The name of the amino acid.
///
#[pyclass]
#[derive(Debug)]
pub struct AminoAcid(rustyms::AminoAcid);

#[pymethods]
impl AminoAcid {
    #[new]
    fn new(name: &str) -> PyResult<Self> {
        match rustyms::AminoAcid::try_from(name) {
            Ok(aa) => Ok(AminoAcid(aa)),
            Err(_) => Err(PyValueError::new_err("Invalid amino acid")),
        }
    }

    fn __str__(&self) -> String {
        self.0.char().to_string()
    }

    fn __repr__(&self) -> String {
        self.to_string()
    }

    /// Molecular formula(s) of the amino acid.
    ///
    /// Returns a list of molecular formulas that are possible for the amino acid symbol.
    ///
    /// Returns
    /// -------
    /// List[MolecularFormula]
    ///
    fn formulas(&self) -> Vec<MolecularFormula> {
        self.0
            .formulas()
            .iter()
            .map(|f| MolecularFormula(f.clone()))
            .collect()
    }

    /// Molecular formula of the amino acid.
    ///
    /// Returns the molecular formula of the amino acid (the first if multiple are possible).
    ///
    /// Returns
    /// -------
    /// List[MolecularFormula]
    ///
    fn formula(&self) -> MolecularFormula {
        MolecularFormula(self.0.formulas().first().unwrap().clone())
    }

    /// Monoisotopic mass(es) of the amino acid.
    ///
    /// Returns a list of monoisotopic masses that are possible for the amino acid symbol.
    ///
    /// Returns
    /// -------
    /// List[float]
    ///
    fn monoisotopic_masses(&self) -> Vec<f64> {
        self.0
            .formulas()
            .iter()
            .map(|f| f.monoisotopic_mass().value)
            .collect()
    }

    /// Monoisotopic mass of the amino acid.
    ///
    /// Returns the monoisotopic mass of the amino acid (the first if multiple are possible).
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn monoisotopic_mass(&self) -> f64 {
        self.0.formulas().first().unwrap().monoisotopic_mass().value
    }
}

impl std::fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Amino acid modification.
///
/// Parameters
/// ----------
/// name : str
///   The name of the modification.
///
#[pyclass]
#[derive(Debug, Clone)]
pub struct Modification(rustyms::Modification);

#[pymethods]
impl Modification {
    #[new]
    fn new(name: &str) -> PyResult<Self> {
        match rustyms::Modification::try_from(name, 0..name.len(), &mut vec![]) {
            Ok(modification) => Ok(Modification(modification.defined().unwrap())),
            Err(_) => Err(PyValueError::new_err("Invalid modification")),
        }
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("Modification('{}')", self.0)
    }

    /// Molecular formula of the modification.
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    fn formula(&self) -> MolecularFormula {
        MolecularFormula(self.0.formula())
    }

    /// Monoisotopic mass of the modification.
    ///
    /// Returns
    /// -------
    /// float
    ///
    fn monoisotopic_mass(&self) -> f64 {
        self.0.formula().monoisotopic_mass().value
    }
}

/// Modification with ambiguous localisation.
///
/// Parameters
/// ----------
/// id : int
///     The id to compare be able to find the other locations where this modifications can be placed.
/// modification : Modification
///     The modification itself.
/// localisation_score : float | None
///     If present the localisation score, meaning the chance/ratio for this modification to show up on this exact spot.
/// group : tuple[str, bool] | None
///     If this is a named group contain the name and track if this is the preferred location or not.
///
#[pyclass]
#[derive(Debug)]
pub struct AmbiguousModification(rustyms::modification::AmbiguousModification);

#[pymethods]
impl AmbiguousModification {
    #[new]
    fn new(
        id: usize,
        modification: Modification,
        localisation_score: Option<f64>,
        group: Option<(String, bool)>,
    ) -> Self {
        AmbiguousModification(rustyms::modification::AmbiguousModification {
            id,
            modification: modification.0,
            localisation_score: localisation_score.map(OrderedFloat),
            group,
        })
    }

    fn __repr__(&self) -> String {
        format!(
            "AmbiguousModification(id={}, modification={}, localisation_score={}, group={:?})",
            self.0.id,
            self.0.modification,
            self.0.localisation_score.unwrap_or_default(),
            match &self.0.group {
                Some((s, b)) => format!("({:?}, {:?})", s, b),
                None => String::new(),
            }
        )
    }

    /// The id to compare be able to find the other locations where this modifications can be placed.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn id(&self) -> usize {
        self.0.id
    }

    /// The modification itself.
    ///
    /// Returns
    /// -------
    /// Modification
    ///
    #[getter]
    fn modification(&self) -> Modification {
        Modification(self.0.modification.clone())
    }

    /// If present the localisation score, meaning the chance/ratio for this modification to show up on this exact spot.
    ///
    /// Returns
    /// -------
    /// float | None
    ///
    #[getter]
    fn localisation_score(&self) -> Option<f64> {
        self.0.localisation_score.map(|x| x.into_inner())
    }

    /// If this is a named group contain the name and track if this is the preferred location or not.
    ///
    /// Returns
    /// -------
    /// tuple[str, bool] | None
    ///
    #[getter]
    fn group(&self) -> Option<(String, bool)> {
        self.0.group.as_ref().map(|(s, b)| (s.to_string(), *b))
    }
}

/// A theoretical fragment of a peptide.
#[pyclass]
#[derive(Debug)]
pub struct Fragment(rustyms::Fragment);

#[pymethods]
impl Fragment {
    fn __repr__(&self) -> String {
        format!(
            "Fragment(formula='{:?}', charge={}, ion='{}', peptide_index={}, neutral_loss='{:?}', label='{}')",
            self.formula(),
            self.charge(),
            self.ion(),
            self.peptide_index(),
            self.neutral_loss(),
            self.label()
        )
    }

    /// The theoretical composition.
    ///
    /// Returns
    /// -------
    /// MolecularFormula
    ///
    #[getter]
    fn formula(&self) -> MolecularFormula {
        MolecularFormula(self.0.formula.clone())
    }

    /// The charge.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn charge(&self) -> i16 {
        self.0.charge.value as i16
    }

    /// All possible annotations for this fragment.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn ion(&self) -> String {
        self.0.ion.to_string() // TODO: Return as exposed Fragment type
    }

    /// The peptide this fragment comes from, saved as the index into the list of peptides in the overarching crate::ComplexPeptide struct.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn peptide_index(&self) -> usize {
        self.0.peptide_index
    }

    /// Any neutral losses applied.
    ///
    /// Returns
    /// -------
    /// str | None
    ///
    #[getter]
    fn neutral_loss(&self) -> Option<String> {
        self.0.neutral_loss.as_ref().map(|nl| nl.to_string())
    }

    /// Additional description for humans.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn label(&self) -> String {
        self.0.label.clone()
    }
}

/// One block in a sequence meaning an amino acid and its accompanying modifications.
#[pyclass]
pub struct SequenceElement(rustyms::SequenceElement);

#[pymethods]
impl SequenceElement {
    fn __repr__(&self) -> String {
        format!("SequenceElement(amino_acid='{}', modifications='{:?}', possible_modifications='{:?}', ambiguous='{:?}')", self.aminoacid(), self.modifications(), self.possible_modifications(), self.ambiguous())
    }

    /// The amino acid.
    ///
    /// Returns
    /// -------
    /// AminoAcid
    ///
    #[getter]
    fn aminoacid(&self) -> AminoAcid {
        AminoAcid(self.0.aminoacid)
    }

    /// All present modifications.
    ///
    /// Returns
    /// -------
    /// list[Modification]
    ///
    #[getter]
    fn modifications(&self) -> Vec<Modification> {
        self.0
            .modifications
            .iter()
            .map(|m| Modification(m.clone()))
            .collect()
    }

    /// All ambiguous modifications (could be placed here or on another position).
    ///
    /// Returns
    /// -------
    /// list[AmbiguousModification]
    ///
    #[getter]
    fn possible_modifications(&self) -> Vec<AmbiguousModification> {
        self.0
            .possible_modifications
            .iter()
            .map(|m| AmbiguousModification(m.clone()))
            .collect()
    }

    /// If this amino acid is part of an ambiguous sequence group `(QA)?` in ProForma
    ///
    /// Returns
    /// -------
    /// int | None
    ///
    #[getter]
    fn ambiguous(&self) -> Option<usize> {
        self.0.ambiguous
    }

    /// Get the molecular formulas for this position with the selected ambiguous modifications, without any global isotype modifications
    ///
    /// Parameters
    /// ----------
    /// selected_ambiguous : int
    ///
    /// Returns
    /// -------
    /// List[MolecularFormula]
    ///
    fn formulas(&self, selected_ambiguous: usize) -> Vec<MolecularFormula> {
        self.0
            .formulas(&[selected_ambiguous])
            .iter()
            .map(|f| MolecularFormula(f.clone()))
            .collect()
    }

    /// Get the molecular formulas for this position with the ambiguous modifications placed on the very first placed (and updating this in `placed`), without any global isotype modifications
    ///
    /// Parameters
    /// ----------
    /// placed : bool
    ///
    /// Returns
    /// -------
    /// List[MolecularFormula]
    ///
    fn formula_greedy(&self, placed: bool) -> Vec<MolecularFormula> {
        self.0
            .formulas_greedy(&mut [placed])
            .iter()
            .map(|f| MolecularFormula(f.clone()))
            .collect()
    }

    /// Get the molecular formulas for this position with all ambiguous modifications, without any global isotype modifications
    ///
    /// Returns
    /// -------
    /// List[MolecularFormula]
    ///
    fn formula_all(&self) -> Vec<MolecularFormula> {
        self.0
            .formulas_all()
            .iter()
            .map(|f| MolecularFormula(f.clone()))
            .collect()
    }
}

/// Fragmentation model enum.
#[pyclass]
enum FragmentationModel {
    All,
    CidHcd,
    Etd,
    Ethcd,
}

/// Helper function to match a [`FragmentationModel`] to a rustyms Model.
fn match_model(model: &FragmentationModel) -> PyResult<rustyms::Model> {
    match model {
        FragmentationModel::All => Ok(rustyms::Model::all()),
        FragmentationModel::CidHcd => Ok(rustyms::Model::cid_hcd()),
        FragmentationModel::Etd => Ok(rustyms::Model::etd()),
        FragmentationModel::Ethcd => Ok(rustyms::Model::ethcd()),
    }
}

/// A peptide with all data as provided by ProForma 2.0.
///
/// Parameters
/// ----------
/// proforma : str
///     The ProForma string.
///
#[pyclass]
#[derive(Clone)]
pub struct LinearPeptide(rustyms::LinearPeptide);

#[pymethods]
impl LinearPeptide {
    /// Create a new peptide from a ProForma string.
    #[new]
    fn new(proforma: &str) -> Self {
        LinearPeptide(rustyms::LinearPeptide::pro_forma(proforma).unwrap())
    }

    fn __str__(&self) -> String {
        self.0.to_string()
    }

    fn __repr__(&self) -> String {
        format!("LinearPeptide({})", self.0)
    }

    fn __len__(&self) -> usize {
        self.0.sequence.len()
    }

    /// Labile modifications, which will not be found in the actual spectrum.
    ///
    /// Returns
    /// -------
    /// list[Modification]
    ///
    #[getter]
    fn labile(&self) -> Vec<Modification> {
        self.0
            .labile
            .iter()
            .map(|x| Modification(x.clone()))
            .collect()
    }

    /// N-terminal modification.
    ///
    /// Returns
    /// -------
    /// Modification | None
    ///
    #[getter]
    fn n_term(&self) -> Option<Modification> {
        self.0.n_term.as_ref().map(|m| Modification(m.clone()))
    }

    /// C-terminal modification.
    ///
    /// Returns
    /// -------
    /// Modification | None
    ///
    #[getter]
    fn c_term(&self) -> Option<Modification> {
        self.0.c_term.as_ref().map(|m| Modification(m.clone()))
    }

    /// Sequence of the peptide including modifications.
    ///
    /// Returns
    /// -------
    /// list[SequenceElement]
    ///
    #[getter]
    fn sequence(&self) -> Vec<SequenceElement> {
        self.0
            .sequence
            .iter()
            .map(|x| SequenceElement(x.clone()))
            .collect()
    }

    /// For each ambiguous modification list all possible positions it can be placed on. Indexed by the ambiguous modification id.
    ///
    /// Returns
    /// -------
    /// list[list[int]]
    ///
    #[getter]
    fn ambiguous_modifications(&self) -> Vec<Vec<usize>> {
        self.0.ambiguous_modifications.clone()
    }

    /// Stripped sequence, meaning the sequence without any modifications.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn stripped_sequence(&self) -> String {
        self.0.sequence.iter().map(|x| x.aminoacid.char()).collect()
    }

    /// The precursor charge of the peptide.
    ///
    /// Returns
    /// -------
    /// int | None
    ///
    #[getter]
    fn charge(&self) -> Option<i16> {
        self.0.charge_carriers.clone().map(|c| c.formula().charge())
    }

    // TODO: Implement when MolecularCharge is exposed upstream.
    // /// The adduct ions, if specified.
    // #[getter]
    // fn charge_carriers(&self) -> Vec<MolecularCharge> {
    //     self.0
    //         .charge_carriers
    //         .iter()
    //         .map(|c| MolecularCharge(c.clone()))
    //         .collect()
    // }

    /// Get a copy of the peptide with its sequence reversed.
    ///
    /// Returns
    /// -------
    /// LinearPeptide
    ///
    fn reverse(&self) -> LinearPeptide {
        LinearPeptide(self.0.reverse())
    }

    /// Gives the formulas for the whole peptide. With the global isotope modifications applied. (Any B/Z will result in multiple possible formulas.)
    ///
    /// Returns
    /// -------
    /// List[MolecularFormula]
    ///
    fn formula(&self) -> Vec<MolecularFormula> {
        self.0
            .formulas()
            .iter()
            .map(|f| MolecularFormula(f.clone()))
            .collect()
    }

    /// Generate the theoretical fragments for this peptide, with the given maximal charge of the fragments, and the given model. With the global isotope modifications applied.
    ///
    /// Parameters
    /// ----------
    /// max_charge : int
    ///     The maximal charge of the fragments.
    /// model : FragmentationModel
    ///     The model to use for the fragmentation.
    ///
    /// Returns
    /// -------
    /// list[Fragment]
    ///   The theoretical fragments.
    ///
    fn generate_theoretical_fragments(
        &self,
        max_charge: i16,
        model: &FragmentationModel,
    ) -> PyResult<Vec<Fragment>> {
        Ok(self
            .0
            .generate_theoretical_fragments(
                rustyms::system::Charge::new::<rustyms::system::e>(max_charge as f64),
                &match_model(model)?,
            )
            .iter()
            .map(|f| Fragment(f.clone()))
            .collect())
    }
}

#[pyclass]
struct RawPeak(rustyms::spectrum::RawPeak);

#[pymethods]
impl RawPeak {
    fn __repr__(&self) -> String {
        format!("RawPeak(mz={}, intensity={})", self.mz(), self.intensity())
    }

    /// The m/z value of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn mz(&self) -> f64 {
        self.0.mz.value
    }

    /// The intensity of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn intensity(&self) -> f64 {
        self.0.intensity.into_inner()
    }
}

impl std::fmt::Display for RawPeak {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "({}, {})", self.mz(), self.intensity())
    }
}

#[pyclass]
/// Represents an annotated peak in a mass spectrometry spectrum.
struct AnnotatedPeak(rustyms::spectrum::AnnotatedPeak);

#[pymethods]
impl AnnotatedPeak {
    fn __repr__(&self) -> String {
        format!(
            "AnnotatedPeak(experimental_mz={}, intensity={}, annotation=[{:?}])",
            self.experimental_mz(),
            self.intensity(),
            self.annotation(),
        )
    }

    /// The experimental m/z value of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn experimental_mz(&self) -> f64 {
        self.0.experimental_mz.value
    }

    /// The intensity of the peak.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn intensity(&self) -> f64 {
        self.0.intensity.into_inner()
    }

    /// All annotations of the peak. Can be empty.
    ///
    /// Returns
    /// -------
    /// list[Fragment]
    ///
    #[getter]
    fn annotation(&self) -> Vec<Fragment> {
        self.0
            .annotation
            .iter()
            .map(|x| Fragment(x.clone()))
            .collect()
    }
}

impl std::fmt::Display for AnnotatedPeak {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "({}, {}, {})",
            self.experimental_mz(),
            self.intensity(),
            self.annotation()
                .iter()
                .map(|x| x.__repr__())
                .collect::<Vec<_>>()
                .join(", ")
        )
    }
}

/// A raw spectrum (meaning not annotated yet)
///
/// Parameters
/// ----------
/// title : str
///     The title of the spectrum.
/// num_scans : int
///     The number of scans.
/// rt : float
///     The retention time.
/// precursor_charge : float
///     The found precursor charge.
/// precursor_mass : float
///     The found precursor mass.
/// mz_array : list[float]
///     The m/z values of the peaks.
/// intensity_array : list[float]
///     The intensities of the peaks.
///
/// Returns
/// -------
/// RawSpectrum
///
#[pyclass]
pub struct RawSpectrum(rustyms::RawSpectrum);

#[pymethods]
impl RawSpectrum {
    /// Create a new raw spectrum.
    ///
    /// Parameters
    /// ----------
    /// title : str
    ///     The title of the spectrum.
    /// num_scans : int
    ///     The number of scans.
    /// rt : float
    ///     The retention time.
    /// precursor_charge : float
    ///     The found precursor charge.
    /// precursor_mass : float
    ///     The found precursor mass.
    /// mz_array : list[float]
    ///     The m/z values of the peaks.
    /// intensity_array : list[float]
    ///     The intensities of the peaks.
    ///
    /// Returns
    /// -------
    /// RawSpectrum
    ///
    #[new]
    fn new(
        title: &str,
        num_scans: u64,
        rt: f64,
        precursor_charge: f64,
        precursor_mass: f64,
        mz_array: Vec<f64>,
        intensity_array: Vec<f64>,
    ) -> Self {
        let mut spectrum = rustyms::RawSpectrum::default();
        spectrum.title = title.to_string();
        spectrum.num_scans = num_scans;
        spectrum.rt = rustyms::system::Time::new::<rustyms::system::s>(rt);
        spectrum.charge = rustyms::system::Charge::new::<rustyms::system::e>(precursor_charge);
        spectrum.mass = rustyms::system::Mass::new::<rustyms::system::dalton>(precursor_mass);

        let peaks = mz_array
            .into_iter()
            .zip(intensity_array)
            .map(|(mz, i)| rustyms::spectrum::RawPeak {
                charge: rustyms::system::Charge::new::<rustyms::system::e>(1.0),
                mz: rustyms::system::MassOverCharge::new::<rustyms::system::mz>(mz),
                intensity: OrderedFloat(i),
            })
            .collect::<Vec<_>>();

        spectrum.extend(peaks);
        RawSpectrum(spectrum)
    }

    fn __repr__(&self) -> String {
        format!(
            "RawSpectrum(title='{}', num_scans={}, rt={}, charge={}, mass={}, spectrum=[{}])",
            self.title(),
            self.num_scans(),
            self.rt(),
            self.charge(),
            self.mass(),
            self.spectrum()
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    /// The title.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn title(&self) -> String {
        self.0.title.clone()
    }

    /// The number of scans.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn num_scans(&self) -> u64 {
        self.0.num_scans
    }

    /// The retention time.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn rt(&self) -> f64 {
        self.0.rt.value
    }

    /// The found precursor charge.
    ///
    /// Returns
    /// -------
    /// float
    #[getter]
    fn charge(&self) -> f64 {
        self.0.charge.value
    }

    /// The found precursor mass.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn mass(&self) -> f64 {
        self.0.mass.value
    }

    /// The peaks of which this spectrum consists.
    ///
    /// Returns
    /// -------
    /// list[RawPeak]
    ///
    #[getter]
    fn spectrum(&self) -> Vec<RawPeak> {
        self.0.clone().into_iter().map(RawPeak).collect()
    }

    /// Annotate this spectrum with the given peptide
    ///
    /// Parameters
    /// ----------
    /// peptide : LinearPeptide
    ///     The peptide to annotate the spectrum with.
    /// model : FragmentationModel
    ///     The model to use for the fragmentation.
    /// mode : MassMode
    ///    The mode to use for the mass.
    ///
    /// Returns
    /// -------
    /// AnnotatedSpectrum
    ///     The annotated spectrum.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the model is not one of the valid models.
    ///
    #[pyo3(signature = (peptide, model, mode=&MassMode::Monoisotopic))]
    fn annotate(
        &self,
        peptide: LinearPeptide,
        model: &FragmentationModel,
        mode: &MassMode,
    ) -> PyResult<AnnotatedSpectrum> {
        let rusty_model = match_model(model)?;
        let fragments = peptide
            .0
            .generate_theoretical_fragments(self.0.charge, &rusty_model);
        Ok(AnnotatedSpectrum(self.0.annotate(
            rustyms::ComplexPeptide::from(peptide.0),
            &fragments,
            &rusty_model,
            match mode {
                MassMode::Monoisotopic => rustyms::MassMode::Monoisotopic,
                MassMode::Average => rustyms::MassMode::Average,
            },
        )))
    }
}

/// An annotated spectrum.
#[pyclass]
pub struct AnnotatedSpectrum(rustyms::AnnotatedSpectrum);

#[pymethods]
impl AnnotatedSpectrum {
    fn __repr__(&self) -> String {
        format!(
            "AnnotatedSpectrum(title='{}', num_scans={}, rt={}, charge={}, mass={}, spectrum=[{}])",
            self.title(),
            self.num_scans(),
            self.rt(),
            self.charge(),
            self.mass(),
            self.spectrum()
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(", ")
        )
    }

    /// The title.
    ///
    /// Returns
    /// -------
    /// str
    ///
    #[getter]
    fn title(&self) -> String {
        self.0.title.clone()
    }

    /// The number of scans.
    ///
    /// Returns
    /// -------
    /// int
    ///
    #[getter]
    fn num_scans(&self) -> u64 {
        self.0.num_scans
    }

    /// The retention time.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn rt(&self) -> f64 {
        self.0.rt.value
    }

    /// The found precursor charge.
    ///
    /// Returns
    /// -------
    /// float
    #[getter]
    fn charge(&self) -> f64 {
        self.0.charge.value
    }

    /// The found precursor mass.
    ///
    /// Returns
    /// -------
    /// float
    ///
    #[getter]
    fn mass(&self) -> f64 {
        self.0.mass.value
    }

    /// The peaks of which this spectrum consists.
    ///
    /// Returns
    /// -------
    /// list[AnnotatedPeak]
    ///
    #[getter]
    fn spectrum(&self) -> Vec<AnnotatedPeak> {
        self.0.clone().into_iter().map(AnnotatedPeak).collect()
    }
}

/// Python bindings to the rustyms library.
#[pymodule]
#[pyo3(name = "rustyms")]
fn rustyms_py03(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<MassMode>()?;
    m.add_class::<Element>()?;
    m.add_class::<MolecularFormula>()?;
    m.add_class::<MolecularCharge>()?;
    m.add_class::<AminoAcid>()?;
    m.add_class::<Modification>()?;
    m.add_class::<AmbiguousModification>()?;
    m.add_class::<Fragment>()?;
    m.add_class::<SequenceElement>()?;
    m.add_class::<FragmentationModel>()?;
    m.add_class::<LinearPeptide>()?;
    m.add_class::<RawPeak>()?;
    m.add_class::<AnnotatedPeak>()?;
    m.add_class::<RawSpectrum>()?;
    m.add_class::<AnnotatedSpectrum>()?;
    Ok(())
}
