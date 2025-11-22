



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



