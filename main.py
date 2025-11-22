
def split_before_each_uppercases(formula):
    parts = []
    current = ""

    for ch in formula:
      
        if ch.isupper() and current:
            parts.append(current)
            current = ch
        else:
            current += ch

    
    if current:
        parts.append(current)

    return parts
 


def split_at_first_digit(formula):
    
    digit_location = 1

    
    for ch in formula[1:]:
        if ch.isdigit():
            break
        digit_location += 1

    
    if digit_location == len(formula):
        return formula, 1

    prefix = formula[:digit_location]
    number = int(formula[digit_location:])
    return prefix, number

    




def count_atoms_in_molecule(molecular_formula):
    """Takes a molecular formula (string) and returns a dictionary of atom counts."""

    parts = split_before_each_uppercases(molecular_formula)
    atom_counts = {}

    for part in parts:
        element, count = split_at_first_digit(part)

        if element not in atom_counts:
            atom_counts[element] = 0

        atom_counts[element] += count

    return atom_counts


def generate_equation_for_element(compounds, coefficients, element):
    """Generates a symbolic equation for the given element from compounds and coefficients.  
    Example: For H in reactants [{'H': 2}, {'O': 4, 'H': 1}], coefficients [a0, a1], returns 2*a0 + a1."""
    equation = 0
    for i, compound in enumerate(compounds):
        if element in compound:
            equation += coefficients[i] * compound[element]
    return equation


def build_equations(reactant_atoms, product_atoms):
    """Builds a list of symbolic equations for each element to balance a chemical reaction.  
    Example: For H2 + O2 -> H2O, returns equations [2*a0 - 2*b0, a1 - b0]."""
    ## coefficients ##
    reactant_coefficients = list(symbols(f'a0:{len(reactant_atoms)}'))
    product_coefficients = list(symbols(f'b0:{len(product_atoms)}')) 
    product_coefficients = product_coefficients[:-1] + [1] # Ensure the last coefficient is 1

    ## equations ##
    equations = []
    for element in ELEMENTS:
        lhs = generate_equation_for_element(reactant_atoms, reactant_coefficients, element)
        rhs = generate_equation_for_element(product_atoms, product_coefficients, element)
        if lhs != 0 or rhs != 0:
            equations.append(Eq(lhs, rhs))

    return equations, reactant_coefficients + product_coefficients[:-1]

from sympy import symbols, Eq, solve as sympy_solve
def my_solve(equations, coefficients):
    """Solves the system of equations for the coefficients of the reaction.  
    Example: For equations [2*a0 - 2*b0, a1 - b0], returns [1.0, 1.0]."""
    solution = sympy_solve(equations, coefficients)

    if len(solution) == len(coefficients):
        coefficient_values = list()
        for coefficient in coefficients:
            coefficient_values.append(float(solution[coefficient]))
        return coefficient_valuesdef parse_chemical_reaction(reaction_equation):
    """Takes a reaction equation (string) and returns reactants and products as lists.  
    Example: 'H2 + O2 -> H2O' → (['H2', 'O2'], ['H2O'])"""
    reaction_equation = reaction_equation.replace(" ", "")  # Remove spaces for easier parsing
    reactants, products = reaction_equation.split("->")
    return reactants.split("+"), products.split("+")

def count_atoms_in_reaction(molecules_list):
    """Takes a list of molecular formulas and returns a list of atom count dictionaries.  
    Example: ['H2', 'O2'] → [{'H': 2}, {'O': 2}]"""
    molecules_atoms_count = []
    for molecule in molecules_list:
        molecules_atoms_count.append(count_atoms_in_molecule(molecule))
    return molecules_atoms_count


def balance_reaction(reaction): #"Fe2O3 + H2 -> Fe + H2O"

    # 1.parse reaction
    reactants, products = parse_chemical_reaction(reaction) # [""Fe2O3", "H2"], ["Fe", "H2O""]
    reactant_atoms = count_atoms_in_reaction(reactants) # [{"Fe":2, "O":1}, {"H":2}]
    product_atoms = count_atoms_in_reaction(products)

    # 2.build equation and solve
    equations, coefficients = build_equations(reactant_atoms, product_atoms)
    coefficients = my_solve(equations, coefficients) + [1]

    return coefficients # [1/3, 1, 2/3, 1]

