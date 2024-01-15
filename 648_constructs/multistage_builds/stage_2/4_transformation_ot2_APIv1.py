from opentrons import labware, instruments, modules, robot
import numpy as np


spotting_tuples=[(('A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1'), ('A1', 'B1', 'C1', 'D1', 'E1', 'F1', 'G1', 'H1'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2'), ('A2', 'B2', 'C2', 'D2', 'E2', 'F2', 'G2', 'H2'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3'), ('A3', 'B3', 'C3', 'D3', 'E3', 'F3', 'G3', 'H3'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4'), ('A4', 'B4', 'C4', 'D4', 'E4', 'F4', 'G4', 'H4'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5'), ('A5', 'B5', 'C5', 'D5', 'E5', 'F5', 'G5', 'H5'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6'), ('A6', 'B6', 'C6', 'D6', 'E6', 'F6', 'G6', 'H6'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7'), ('A7', 'B7', 'C7', 'D7', 'E7', 'F7', 'G7', 'H7'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8'), ('A8', 'B8', 'C8', 'D8', 'E8', 'F8', 'G8', 'H8'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9'), ('A9', 'B9', 'C9', 'D9', 'E9', 'F9', 'G9', 'H9'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10'), ('A10', 'B10', 'C10', 'D10', 'E10', 'F10', 'G10', 'H10'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11'), ('A11', 'B11', 'C11', 'D11', 'E11', 'F11', 'G11', 'H11'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12'), ('A12', 'B12', 'C12', 'D12', 'E12', 'F12', 'G12', 'H12'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A13', 'B13', 'C13', 'D13', 'E13', 'F13', 'G13', 'H13'), ('A13', 'B13', 'C13', 'D13', 'E13', 'F13', 'G13', 'H13'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A14', 'B14', 'C14', 'D14', 'E14', 'F14', 'G14', 'H14'), ('A14', 'B14', 'C14', 'D14', 'E14', 'F14', 'G14', 'H14'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A15', 'B15', 'C15', 'D15', 'E15', 'F15', 'G15', 'H15'), ('A15', 'B15', 'C15', 'D15', 'E15', 'F15', 'G15', 'H15'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A16', 'B16', 'C16', 'D16', 'E16', 'F16', 'G16', 'H16'), ('A16', 'B16', 'C16', 'D16', 'E16', 'F16', 'G16', 'H16'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A17', 'B17', 'C17', 'D17', 'E17', 'F17', 'G17', 'H17'), ('A17', 'B17', 'C17', 'D17', 'E17', 'F17', 'G17', 'H17'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A18', 'B18', 'C18', 'D18', 'E18', 'F18', 'G18', 'H18'), ('A18', 'B18', 'C18', 'D18', 'E18', 'F18', 'G18', 'H18'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A19', 'B19', 'C19', 'D19', 'E19', 'F19', 'G19', 'H19'), ('A19', 'B19', 'C19', 'D19', 'E19', 'F19', 'G19', 'H19'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A20', 'B20', 'C20', 'D20', 'E20', 'F20', 'G20', 'H20'), ('A20', 'B20', 'C20', 'D20', 'E20', 'F20', 'G20', 'H20'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A21', 'B21', 'C21', 'D21', 'E21', 'F21', 'G21', 'H21'), ('A21', 'B21', 'C21', 'D21', 'E21', 'F21', 'G21', 'H21'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A22', 'B22', 'C22', 'D22', 'E22', 'F22', 'G22', 'H22'), ('A22', 'B22', 'C22', 'D22', 'E22', 'F22', 'G22', 'H22'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A23', 'B23', 'C23', 'D23', 'E23', 'F23', 'G23', 'H23'), ('A23', 'B23', 'C23', 'D23', 'E23', 'F23', 'G23', 'H23'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A24', 'B24', 'C24', 'D24', 'E24', 'F24', 'G24', 'H24'), ('A24', 'B24', 'C24', 'D24', 'E24', 'F24', 'G24', 'H24'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A25', 'B25', 'C25', 'D25', 'E25', 'F25', 'G25', 'H25'), ('A25', 'B25', 'C25', 'D25', 'E25', 'F25', 'G25', 'H25'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A26', 'B26', 'C26', 'D26', 'E26', 'F26', 'G26', 'H26'), ('A26', 'B26', 'C26', 'D26', 'E26', 'F26', 'G26', 'H26'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A27', 'B27', 'C27', 'D27', 'E27', 'F27', 'G27', 'H27'), ('A27', 'B27', 'C27', 'D27', 'E27', 'F27', 'G27', 'H27'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A28', 'B28', 'C28', 'D28', 'E28', 'F28', 'G28', 'H28'), ('A28', 'B28', 'C28', 'D28', 'E28', 'F28', 'G28', 'H28'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A29', 'B29', 'C29', 'D29', 'E29', 'F29', 'G29', 'H29'), ('A29', 'B29', 'C29', 'D29', 'E29', 'F29', 'G29', 'H29'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A30', 'B30', 'C30', 'D30', 'E30', 'F30', 'G30', 'H30'), ('A30', 'B30', 'C30', 'D30', 'E30', 'F30', 'G30', 'H30'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A31', 'B31', 'C31', 'D31', 'E31', 'F31', 'G31', 'H31'), ('A31', 'B31', 'C31', 'D31', 'E31', 'F31', 'G31', 'H31'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A32', 'B32', 'C32', 'D32', 'E32', 'F32', 'G32', 'H32'), ('A32', 'B32', 'C32', 'D32', 'E32', 'F32', 'G32', 'H32'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A33', 'B33', 'C33', 'D33', 'E33', 'F33', 'G33', 'H33'), ('A33', 'B33', 'C33', 'D33', 'E33', 'F33', 'G33', 'H33'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A34', 'B34', 'C34', 'D34', 'E34', 'F34', 'G34', 'H34'), ('A34', 'B34', 'C34', 'D34', 'E34', 'F34', 'G34', 'H34'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A35', 'B35', 'C35', 'D35', 'E35', 'F35', 'G35', 'H35'), ('A35', 'B35', 'C35', 'D35', 'E35', 'F35', 'G35', 'H35'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A36', 'B36', 'C36', 'D36', 'E36', 'F36', 'G36', 'H36'), ('A36', 'B36', 'C36', 'D36', 'E36', 'F36', 'G36', 'H36'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A37', 'B37', 'C37', 'D37', 'E37', 'F37', 'G37', 'H37'), ('A37', 'B37', 'C37', 'D37', 'E37', 'F37', 'G37', 'H37'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A38', 'B38', 'C38', 'D38', 'E38', 'F38', 'G38', 'H38'), ('A38', 'B38', 'C38', 'D38', 'E38', 'F38', 'G38', 'H38'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A39', 'B39', 'C39', 'D39', 'E39', 'F39', 'G39', 'H39'), ('A39', 'B39', 'C39', 'D39', 'E39', 'F39', 'G39', 'H39'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A40', 'B40', 'C40', 'D40', 'E40', 'F40', 'G40', 'H40'), ('A40', 'B40', 'C40', 'D40', 'E40', 'F40', 'G40', 'H40'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A41', 'B41', 'C41', 'D41', 'E41', 'F41', 'G41', 'H41'), ('A41', 'B41', 'C41', 'D41', 'E41', 'F41', 'G41', 'H41'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A42', 'B42', 'C42', 'D42', 'E42', 'F42', 'G42', 'H42'), ('A42', 'B42', 'C42', 'D42', 'E42', 'F42', 'G42', 'H42'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A43', 'B43', 'C43', 'D43', 'E43', 'F43', 'G43', 'H43'), ('A43', 'B43', 'C43', 'D43', 'E43', 'F43', 'G43', 'H43'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A44', 'B44', 'C44', 'D44', 'E44', 'F44', 'G44', 'H44'), ('A44', 'B44', 'C44', 'D44', 'E44', 'F44', 'G44', 'H44'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A45', 'B45', 'C45', 'D45', 'E45', 'F45', 'G45', 'H45'), ('A45', 'B45', 'C45', 'D45', 'E45', 'F45', 'G45', 'H45'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A46', 'B46', 'C46', 'D46', 'E46', 'F46', 'G46', 'H46'), ('A46', 'B46', 'C46', 'D46', 'E46', 'F46', 'G46', 'H46'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A47', 'B47', 'C47', 'D47', 'E47', 'F47', 'G47', 'H47'), ('A47', 'B47', 'C47', 'D47', 'E47', 'F47', 'G47', 'H47'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A48', 'B48', 'C48', 'D48', 'E48', 'F48', 'G48', 'H48'), ('A48', 'B48', 'C48', 'D48', 'E48', 'F48', 'G48', 'H48'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A49', 'B49', 'C49', 'D49', 'E49', 'F49', 'G49', 'H49'), ('A49', 'B49', 'C49', 'D49', 'E49', 'F49', 'G49', 'H49'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A50', 'B50', 'C50', 'D50', 'E50', 'F50', 'G50', 'H50'), ('A50', 'B50', 'C50', 'D50', 'E50', 'F50', 'G50', 'H50'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A51', 'B51', 'C51', 'D51', 'E51', 'F51', 'G51', 'H51'), ('A51', 'B51', 'C51', 'D51', 'E51', 'F51', 'G51', 'H51'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A52', 'B52', 'C52', 'D52', 'E52', 'F52', 'G52', 'H52'), ('A52', 'B52', 'C52', 'D52', 'E52', 'F52', 'G52', 'H52'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A53', 'B53', 'C53', 'D53', 'E53', 'F53', 'G53', 'H53'), ('A53', 'B53', 'C53', 'D53', 'E53', 'F53', 'G53', 'H53'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A54', 'B54', 'C54', 'D54', 'E54', 'F54', 'G54', 'H54'), ('A54', 'B54', 'C54', 'D54', 'E54', 'F54', 'G54', 'H54'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A55', 'B55', 'C55', 'D55', 'E55', 'F55', 'G55', 'H55'), ('A55', 'B55', 'C55', 'D55', 'E55', 'F55', 'G55', 'H55'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A56', 'B56', 'C56', 'D56', 'E56', 'F56', 'G56', 'H56'), ('A56', 'B56', 'C56', 'D56', 'E56', 'F56', 'G56', 'H56'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A57', 'B57', 'C57', 'D57', 'E57', 'F57', 'G57', 'H57'), ('A57', 'B57', 'C57', 'D57', 'E57', 'F57', 'G57', 'H57'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A58', 'B58', 'C58', 'D58', 'E58', 'F58', 'G58', 'H58'), ('A58', 'B58', 'C58', 'D58', 'E58', 'F58', 'G58', 'H58'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A59', 'B59', 'C59', 'D59', 'E59', 'F59', 'G59', 'H59'), ('A59', 'B59', 'C59', 'D59', 'E59', 'F59', 'G59', 'H59'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A60', 'B60', 'C60', 'D60', 'E60', 'F60', 'G60', 'H60'), ('A60', 'B60', 'C60', 'D60', 'E60', 'F60', 'G60', 'H60'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A61', 'B61', 'C61', 'D61', 'E61', 'F61', 'G61', 'H61'), ('A61', 'B61', 'C61', 'D61', 'E61', 'F61', 'G61', 'H61'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A62', 'B62', 'C62', 'D62', 'E62', 'F62', 'G62', 'H62'), ('A62', 'B62', 'C62', 'D62', 'E62', 'F62', 'G62', 'H62'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A63', 'B63', 'C63', 'D63', 'E63', 'F63', 'G63', 'H63'), ('A63', 'B63', 'C63', 'D63', 'E63', 'F63', 'G63', 'H63'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A64', 'B64', 'C64', 'D64', 'E64', 'F64', 'G64', 'H64'), ('A64', 'B64', 'C64', 'D64', 'E64', 'F64', 'G64', 'H64'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A65', 'B65', 'C65', 'D65', 'E65', 'F65', 'G65', 'H65'), ('A65', 'B65', 'C65', 'D65', 'E65', 'F65', 'G65', 'H65'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A66', 'B66', 'C66', 'D66', 'E66', 'F66', 'G66', 'H66'), ('A66', 'B66', 'C66', 'D66', 'E66', 'F66', 'G66', 'H66'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A67', 'B67', 'C67', 'D67', 'E67', 'F67', 'G67', 'H67'), ('A67', 'B67', 'C67', 'D67', 'E67', 'F67', 'G67', 'H67'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A68', 'B68', 'C68', 'D68', 'E68', 'F68', 'G68', 'H68'), ('A68', 'B68', 'C68', 'D68', 'E68', 'F68', 'G68', 'H68'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A69', 'B69', 'C69', 'D69', 'E69', 'F69', 'G69', 'H69'), ('A69', 'B69', 'C69', 'D69', 'E69', 'F69', 'G69', 'H69'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A70', 'B70', 'C70', 'D70', 'E70', 'F70', 'G70', 'H70'), ('A70', 'B70', 'C70', 'D70', 'E70', 'F70', 'G70', 'H70'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A71', 'B71', 'C71', 'D71', 'E71', 'F71', 'G71', 'H71'), ('A71', 'B71', 'C71', 'D71', 'E71', 'F71', 'G71', 'H71'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A72', 'B72', 'C72', 'D72', 'E72', 'F72', 'G72', 'H72'), ('A72', 'B72', 'C72', 'D72', 'E72', 'F72', 'G72', 'H72'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A73', 'B73', 'C73', 'D73', 'E73', 'F73', 'G73', 'H73'), ('A73', 'B73', 'C73', 'D73', 'E73', 'F73', 'G73', 'H73'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A74', 'B74', 'C74', 'D74', 'E74', 'F74', 'G74', 'H74'), ('A74', 'B74', 'C74', 'D74', 'E74', 'F74', 'G74', 'H74'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A75', 'B75', 'C75', 'D75', 'E75', 'F75', 'G75', 'H75'), ('A75', 'B75', 'C75', 'D75', 'E75', 'F75', 'G75', 'H75'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A76', 'B76', 'C76', 'D76', 'E76', 'F76', 'G76', 'H76'), ('A76', 'B76', 'C76', 'D76', 'E76', 'F76', 'G76', 'H76'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A77', 'B77', 'C77', 'D77', 'E77', 'F77', 'G77', 'H77'), ('A77', 'B77', 'C77', 'D77', 'E77', 'F77', 'G77', 'H77'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A78', 'B78', 'C78', 'D78', 'E78', 'F78', 'G78', 'H78'), ('A78', 'B78', 'C78', 'D78', 'E78', 'F78', 'G78', 'H78'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A79', 'B79', 'C79', 'D79', 'E79', 'F79', 'G79', 'H79'), ('A79', 'B79', 'C79', 'D79', 'E79', 'F79', 'G79', 'H79'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A80', 'B80', 'C80', 'D80', 'E80', 'F80', 'G80', 'H80'), ('A80', 'B80', 'C80', 'D80', 'E80', 'F80', 'G80', 'H80'), (5, 5, 5, 5, 5, 5, 5, 5)), (('A81', 'B81', 'C81', 'D81', 'E81', 'F81', 'G81', 'H81'), ('A81', 'B81', 'C81', 'D81', 'E81', 'F81', 'G81', 'H81'), (5, 5, 5, 5, 5, 5, 5, 5))]
soc_well='A1'


def generate_transformation_wells(spotting_tuples):
    """Evaluates spotting_tuples and returns transformation wells.

    Args:
    spotting_tuples (list): Sets of spotting reactions are given 
        in the form: ((source wells), (target wells), (spotting volumes)). 
        Each unique transformation well is resuspended once prior to spotting.

    """
    wells = []
    for spotting_tuple in spotting_tuples:
        for source_well in spotting_tuple[0]:
            wells.append(source_well)
    transformation_wells = [well for i, well in enumerate(
        wells) if wells.index(well) == i]
    return transformation_wells


def tiprack_slots(
        spotting_tuples,
        max_spot_vol=5):
    """Calculates p10 and p300 tiprack slots required.

    Args:
    spotting_tuples (list): Sets of spotting reactions are given 
        in the form: ((source wells), (target wells), (spotting volumes)). 
        Each unique transformation well is resuspended once prior to spotting.
    max_spot_vol (float): Maximum volume that is spotted per spot reaction.

    """
    # Reactions' number
    transformation_reactions = len(
        generate_transformation_wells(spotting_tuples))
    spotting_reactions = 0
    for spotting_tuple in spotting_tuples:
        spots = np.array(spotting_tuple[2])/max_spot_vol
        np.ceil(spots)
        spotting_reactions = spotting_reactions + int(np.sum(spots))

    # p10 tiprack slots
    p10_tips = transformation_reactions + spotting_reactions
    p10_tiprack_slots = p10_tips // 96 + 1 if p10_tips % 96 > 0 else p10_tips / 96

    # p300 tiprack slots
    p300_tips = transformation_reactions + spotting_reactions
    p300_tiprack_slots = p300_tips // 96 + \
        1 if p300_tips % 96 > 0 else p300_tips / 96
    return int(p10_tiprack_slots), int(p300_tiprack_slots)


def transformation_setup(transformation_wells):
    """Sets up transformation reactions

    Args:
    transformation_wells (list). 

    """

    # Constants
    TEMP = 4  # Incubation temperature.
    ASSEMBLY_VOL = 5  # Volume of final assembly added to competent cells.
    MIX_SETTINGS = (4, 5)  # Mix after setting during final assembly transfers.
    INCUBATION_TIME = 20  # Cells and final assembly incubation time.

    # Set temperature deck to 4 °C and load competent cells
    tempdeck.set_temperature(TEMP)
    tempdeck.wait_for_temp()
    robot.pause()
    robot.comment('Load competent cells, uncap and resume run')

    # Transfer final assemblies
    p10_pipette.transfer(ASSEMBLY_VOL,
                         assembly_plate.wells(transformation_wells),
                         transformation_plate.wells(
                             transformation_wells), new_tip='always',
                         mix_after=(MIX_SETTINGS))

    # Incubate for 20 minutes and remove competent cells for heat shock
    p10_pipette.delay(minutes=INCUBATION_TIME)
    robot.pause()
    robot.comment(
        'Remove transformation reactions, conduct heatshock and replace.')


def phase_switch(comment='Remove final assembly plate. Introduce agar tray and deep well plate containing SOC media. Resume run.'):
    """Function pauses run enabling addition/removal of labware.

    Args:
    comment (str): string to be displayed during run following pause.

    """
    robot.pause()
    robot.comment(comment)


def outgrowth(
        cols,
        soc_well):
    """Outgrows transformed cells.

    Args:
    cols (list of str): list of cols in transformation plate containing samples.
    soc_well (str): Well containing SOC media in relevant plate.

    """

    # Constants
    SOC_VOL = 125
    SOC_MIX_SETTINGS = (4, 50)
    TEMP = 37
    OUTGROWTH_TIME = 60
    SOC_ASPIRATION_RATE = 25
    P300_DEFAULT_ASPIRATION_RATE = 150

    # Define wells
    transformation_cols = transformation_plate.cols(cols)
    soc = soc_plate.wells(soc_well)

    # Add SOC to transformed cells
    p300_pipette.set_flow_rate(aspirate=SOC_ASPIRATION_RATE)
    p300_pipette.transfer(SOC_VOL, soc, transformation_cols,
                          new_tip='always', mix_after=SOC_MIX_SETTINGS)
    p300_pipette.set_flow_rate(aspirate=P300_DEFAULT_ASPIRATION_RATE)

    # Incubate for 1 hour at 37 °C
    tempdeck.set_temperature(TEMP)
    tempdeck.wait_for_temp()
    p300_pipette.delay(minutes=OUTGROWTH_TIME)
    tempdeck.deactivate()


def spotting_cols(spotting_tuples):
    """Evaluates spotting_tuples and returns unique cols (str) 
    associated with each spotting_tuple's source wells.

    Args:
    spotting_tuples (list): Sets of spotting reactions are given 
        in the form: ((source wells), (target wells), (spotting volumes)). 
        Each unique transformation well is resuspended once prior to spotting.

    """
    cols_list = []
    for spotting_tuple in spotting_tuples:
        source_wells_cols = [source_well[1:]
                             for source_well in spotting_tuple[0]]
        unique_cols = [col for i, col in enumerate(
            source_wells_cols) if source_wells_cols.index(col) == i]
        cols_list.append(unique_cols)
    return cols_list


def spot_transformations(
        spotting_tuples,
        dead_vol=2,
        spotting_dispense_rate=0.025,
        stabbing_depth=10,
        max_spot_vol=5):
    """Spots transformation reactions.

    Args:
    spotting_tuples (list): Sets of spotting reactions are given 
        in the form: ((source wells), (target wells), (spotting volumes)). 
        Each unique source well is resuspended once prior to spotting.
    dead_vol (float): Dead volume aspirated during spotting.
    spotting_dispense_rate (float): Rate p10_pipette dispenses at during spotting.
    stabbing_depth (float): Depth p10_pipette moves into agar during spotting.
    max_spot_vol (float): Maximum volume that is spotted per spot reaction. 

    """

    def spot(
            source,
            target,
            spot_vol):
        """Spots an individual reaction using the p10 pipette.

        Args:
        source (str): Well containing the transformation reaction to be spotted.
        target (str): Well transformation reaction is to be spotted to.
        spot_vol (float): Volume of transformation reaction to be spotted (uL).  

        """
        # Constants
        DEFAULT_HEAD_SPEED = {'x': 400, 'y': 400,
                              'z': 125, 'a': 125}
        SPOT_HEAD_SPEED = {'x': 400, 'y': 400, 'z': 125,
                           'a': 125 // 4}
        DISPENSING_HEIGHT = 5
        SAFE_HEIGHT = 15  # height avoids collision with agar tray.

        # Spot
        p10_pipette.pick_up_tip()
        p10_pipette.aspirate(spot_vol + dead_vol, source)
        p10_pipette.move_to(target.top(SAFE_HEIGHT))
        p10_pipette.move_to(target.top(DISPENSING_HEIGHT))
        p10_pipette.dispense(volume=spot_vol, rate=spotting_dispense_rate)
        robot.head_speed(combined_speed=max(
            SPOT_HEAD_SPEED.values()), **SPOT_HEAD_SPEED)
        p10_pipette.move_to(target.top(-1 * stabbing_depth))
        robot.head_speed(combined_speed=max(
            DEFAULT_HEAD_SPEED.values()), **DEFAULT_HEAD_SPEED)
        p10_pipette.move_to(target.top(SAFE_HEIGHT))

        # Dispose of dead volume and tip
        p10_pipette.dispense(dead_vol, spotting_waste)
        p10_pipette.blow_out()
        p10_pipette.drop_tip()

    def spot_tuple(spotting_tuple):
        """Spots all reactions defined by the spotting tuple. Requires the function spot.

            Args:
            spotting_tuple (tuple): Spotting reactions given in the form: (source wells), (target wells), (spotting volumes).

        """
        source_wells = spotting_tuple[0]
        target_wells = spotting_tuple[1]
        spot_vols = list(spotting_tuple[2])
        while max(spot_vols) > 0:
            for index, spot_vol in enumerate(spot_vols):
                if spot_vol == 0:
                    pass
                else:
                    vol = spot_vol if spot_vol <= max_spot_vol else max_spot_vol
                    spot(transformation_plate.wells(source_wells[index]),
                         agar_plate.wells(target_wells[index]), vol)
                    spot_vols[index] = spot_vols[index] - vol

    # Constants
    TRANSFORMATION_MIX_SETTINGS = [4, 50]

    # Spot transformation reactions
    for spotting_tuple in spotting_tuples:
        source_wells_cols = [source_well[1:]
                             for source_well in spotting_tuple[0]]
        unique_cols = [col for i, col in enumerate(
            source_wells_cols) if source_wells_cols.index(col) == i]
        for col in unique_cols:
            p300_pipette.pick_up_tip()
            p300_pipette.mix(TRANSFORMATION_MIX_SETTINGS[0],
                             TRANSFORMATION_MIX_SETTINGS[1],
                             transformation_plate.cols(col))
            p300_pipette.drop_tip()
        spot_tuple(spotting_tuple)


# Run protocol

# Constants
CANDIDATE_P10_SLOTS = ['9', '2', '5']
CANDIDATE_P300_SLOTS = ['3', '6']
P10_TIPRACK_TYPE = 'tiprack-10ul'
P300_TIPRACK_TYPE = 'opentrons_96_tiprack_300ul'
P10_MOUNT = 'right'
P300_MOUNT = 'left'
ASSEMBLY_PLATE_TYPE = '4ti-0960_FrameStar'
ASSEMBLY_PLATE_SLOT = '8'
TEMPDECK_SLOT = '10'
TRANSFORMATION_PLATE_TYPE = 'Eppendorf_30133366_plate_96'
SOC_PLATE_TYPE = '4ti0136_96_deep-well'
SOC_PLATE_SLOT = '7'
TUBE_RACK_TYPE = 'tube-rack_E1415-1500'
TUBE_RACK_SLOT = '11'
SPOTTING_WASTE_WELL = 'A1'
AGAR_PLATE_TYPE = 'Nunc_Omnitray'
AGAR_PLATE_SLOT = '1'

# Tiprack slots
p10_p300_tiprack_slots = tiprack_slots(spotting_tuples)
p10_slots = CANDIDATE_P10_SLOTS[
    :p10_p300_tiprack_slots[0]]
p300_slots = CANDIDATE_P300_SLOTS[
    :p10_p300_tiprack_slots[1]]

# Define labware
p10_tipracks = [labware.load(P10_TIPRACK_TYPE, slot)
                for slot in p10_slots]
p300_tipracks = [labware.load(P300_TIPRACK_TYPE, slot)
                 for slot in p300_slots]
p10_pipette = instruments.P10_Single(
    mount=P10_MOUNT, tip_racks=p10_tipracks)
p300_pipette = instruments.P300_Multi(
    mount=P300_MOUNT, tip_racks=p300_tipracks)

assembly_plate = labware.load(ASSEMBLY_PLATE_TYPE, ASSEMBLY_PLATE_SLOT)
tempdeck = modules.load('tempdeck', TEMPDECK_SLOT)
transformation_plate = labware.load(TRANSFORMATION_PLATE_TYPE,
                                    TEMPDECK_SLOT, share=True)
soc_plate = labware.load(SOC_PLATE_TYPE, SOC_PLATE_SLOT)
tube_rack = labware.load(TUBE_RACK_TYPE, TUBE_RACK_SLOT)
spotting_waste = tube_rack.wells(SPOTTING_WASTE_WELL)
agar_plate = labware.load(AGAR_PLATE_TYPE, AGAR_PLATE_SLOT)

# Register agar_plate for calibration
p10_pipette.transfer(1, agar_plate.wells(
    'A1'), agar_plate.wells('H12'), trash=False)
p10_pipette.start_at_tip(p10_tipracks[0][0])

# Run functions
transformation_setup(generate_transformation_wells(spotting_tuples))
phase_switch()
spotting_tuples_cols = [col for cols in spotting_cols(
    spotting_tuples) for col in cols]
unique_cols = [col for i, col in enumerate(
    spotting_tuples_cols) if spotting_tuples_cols.index(col) == i]
outgrowth(unique_cols, soc_well=soc_well)
spot_transformations(spotting_tuples)

for c in robot.commands():
    print(c)