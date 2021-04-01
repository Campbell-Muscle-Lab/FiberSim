"""
@author: Dylan Colli

Half-sarcomere classes.
"""

class ThickFilament():
  """Thick filament"""
  def __init__(self, d_t):
    self.id = d_t["thick_id"]  # Change for new format.
    self.m_y = d_t["m_y"]
    self.m_z = d_t["m_z"]
    self.m_no_of_cbs = d_t["m_no_of_cbs"]
    self.m_k_stiff = d_t["m_k_stiff"]
    self.m_inter_crown_rest_length = d_t["m_inter_crown_rest_length"]
    self.m_lambda = d_t["m_lambda"]   
    self.nearest_actin_filaments = d_t["nearest_actin_filaments"]
    self.cb_x = d_t["cb_x"]
    self.cb_angle = d_t["cb_angle"]
    self.cb_state = d_t["cb_state"]
    self.cb_isoform = d_t["cb_iso"]
    self.cb_bound_to_a_f = d_t["cb_bound_to_a_f"]
    self.cb_bound_to_a_n = d_t["cb_bound_to_a_n"]
    self.cb_nearest_a_f = d_t["cb_nearest_a_f"]
    self.cb_nearest_a_n = d_t["cb_nearest_a_n"]
    self.pc_node_index = d_t["pc_node_index"]
    self.pc_angle = d_t["pc_angle"]
    self.pc_state = None
    self.pc_phos = None
    self.pc_bound_to_a_f = None
    self.pc_bound_to_a_n = None
    self.pc_nearest_a_f = None
    self.pc_nearest_a_n = None


class ThinFilament():
  """Thin filament"""
  def __init__(self, d_t):
    self.id = d_t["thin_id"]
    self.a_y = d_t["a_y"]
    self.a_z = d_t["a_z"]
    self.a_no_of_bs = d_t["a_no_of_bs"]
    self.a_k_stiff = d_t["a_k_stiff"]
    self.a_inter_bs_rest_length = d_t["a_inter_bs_rest_length"]
    self.bs_x = d_t["bs_x"]
    self.bs_angle = d_t["bs_angle"]
    self.bs_unit = d_t["bs_angle"]
    self.bs_state = d_t["bs_state"]
    self.bs_isoform = d_t["bs_isoform"]
    self.bound_to_m_f = d_t["bound_to_m_f"]
    self.bound_to_m_n = d_t["bound_to_m_n"]
    self.nearest_m_f = d_t["nearest_m_f"]
    self.nearest_m_n = d_t["nearest_m_n"]


class HalfSarcomere():
  """Half-sarcomere"""
  def __init__(self, dump_dict):
    head = dump_dict["hs_data"]
    self.hs_id = head["hs_id"]
    self.time = head["time"]
    self.hs_length = head["hs_length"]
    self.hs_force = head["hs_force"]
    self.pCa = head["pCa"]
    self.m_nodes_per_thick_filament = head["m_nodes_per_thick_filament"]
    self.a_nodes_per_thin_filament = head["a_nodes_per_thin_filament"]
    self.thick = []
    self.thin = []

    # Construct all of the thick filaments.
    for m_dict in head["thick"]:
      self.thick.append(ThickFilament(m_dict))

    # Construct all of the thin filaments.
    for a_dict in head["thin"]:
      self.thin.append(ThinFilament(a_dict))