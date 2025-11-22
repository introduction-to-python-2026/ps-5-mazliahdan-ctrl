


def split_before_uppercases(formula):
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

def split_at_digit(formula):
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
    """Takes a molecular formula (string) and returns a dictionary of atom counts.  
    Example: 'H2O' → {'H': 2, 'O': 1}"""

    # Step 1: Initialize an empty dictionary to store atom counts

    for atom in split_by_capitals(molecular_formula):
        atom_name, atom_count = split_at_number(atom)
        
        # Step 2: Update the dictionary with the atom name and count

    # Step 3: Return the completed dictionary



def parse_chemical_reaction(reaction_equation):
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
