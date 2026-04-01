# optimised by Sofie and Trisha :))
# Modern UI Enhancement

from PySide6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QLabel, QComboBox, QPushButton,
    QFileDialog, QFormLayout, QScrollArea, QFrame, QHBoxLayout, QLineEdit
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap, QFont
import sys
import os

# ---------------------------------------
# LABWARE + HARDWARE DEFINITIONS
# ---------------------------------------

OT2_SINGLE_PIPETTES = ["p20_single_gen2"]
FLEX_SINGLE_PIPETTES = ["Flex_1channel_50"]

OT2_MULTI_PIPETTES = ["p300_multi_gen2"]
FLEX_MULTI_PIPETTES = ["Flex_8channel_1000"]

OT2_DEFAULTS = {
    "single_pipette": "p20_single_gen2",
    "multi_pipette": "p300_multi_gen2",
    "single_mount": "right",
    "multi_mount": "left",
    "thermocycler": "thermocyclerModuleV2",
    "magdeck": "magnetic module gen2",
}

FLEX_DEFAULTS = {
    "single_pipette": "Flex_1channel_50",
    "multi_pipette": "Flex_8channel_1000",
    "single_mount": "right",
    "multi_mount": "left",
    "thermocycler": "thermocyclerModuleV2",
    "magdeck": "magneticBlockV1",
}

OT2_LABWARE = {
    "24_tuberack_1500ul": ["e14151500starlab_24_tuberack_1500ul"],
    "tiprack_20ul": ["opentrons_96_tiprack_20ul"],
    "tiprack_300ul": ["opentrons_96_tiprack_300ul"],
    "clip_source_plate": ["4ti0960rig_96_wellplate_200ul"],
    "clip_plate": ["4ti0960rig_96_wellplate_200ul"],
    "mix_plate": ["4ti0960rig_96_wellplate_200ul"],
    "mag_plate": ["4ti0960rig_96_wellplate_200ul"],
    "final_assembly_plate": ["4ti0960rig_96_wellplate_200ul"],
    "transform_plate": ["4ti0960rig_96_wellplate_200ul"],
    "agar_plate": ["4ti0960rig_96_wellplate_200ul"],
    "12_reservoir_21000ul": ["4ti0131_12_reservoir_21000ul"],
    "96_deepwellplate_2ml": ["4ti0136_96_wellplate_2200ul"],
    "12_corning_wellplate": ["corning_12_wellplate_6.9ml_flat"]
}

FLEX_LABWARE = {
    "24_tuberack_1500ul": ["e14151500starlab_24_tuberack_1500ul"],
    "tiprack_20ul": ["opentrons_flex_96_tiprack_50ul"],
    "tiprack_300ul": ["opentrons_flex_96_tiprack_200ul"],
    "flex_96_tiprack_200ul": ["opentrons_flex_96_tiprack_200ul"],
    "flex_96_tiprack_1000ul": ["opentrons_flex_96_tiprack_1000ul"],
    "clip_source_plate": ["4ti0960rig_96_wellplate_200ul"],
    "clip_plate": ["4ti0960rig_96_wellplate_200ul"],
    "mix_plate": ["4ti0960rig_96_wellplate_200ul"],
    "mag_plate": ["4ti0960rig_96_wellplate_200ul"],
    "final_assembly_plate": ["4ti0960rig_96_wellplate_200ul"],
    "transform_plate": ["4ti0960rig_96_wellplate_200ul"],
    "agar_plate": ["4ti0960rig_96_wellplate_200ul"],
    "12_reservoir_21000ul": ["4ti0131_12_reservoir_21000ul"],
    "96_deepwellplate_2ml": ["4ti0136_96_wellplate_2200ul"],
    "12_corning_wellplate": ["corning_12_wellplate_6.9ml_flat"]
}

ALL_LABWARE_KEYS = list(dict.fromkeys(list(OT2_LABWARE.keys()) + list(FLEX_LABWARE.keys())))

# ---------------------------------------
# MODERN STYLESHEET
# ---------------------------------------

MODERN_STYLESHEET = """
QWidget {
    background-color: #E8F4F8;
    color: #0A3D62;
    font-family: 'Helvetica Neue', Arial, sans-serif;
    font-size: 13px;
}

QLabel {
    color: #0A3D62;
    padding: 2px;
    background: transparent;
}



QComboBox {
    background-color: #FFFFFF;
    border: 2px solid #5DADE2;
    border-radius: 6px;
    padding: 8px 12px;
    color: #0A3D62;
    min-height: 25px;
    min-width: 200px;
    font-weight: 500;
}

QComboBox:hover {
    border: 2px solid #3498DB;
    background-color: #F8FBFD;
}

QComboBox:focus {
    border: 2px solid #2E86C1;
}

QComboBox::drop-down {
    border: none;
    width: 30px;
}

QComboBox::down-arrow {
    image: none;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 6px solid #3498DB;
    margin-right: 8px;
}


QComboBox QAbstractItemView {
    background-color: #FFFFFF;
    border: 2px solid #5DADE2;
    selection-background-color: #AED6F1;
    selection-color: #0A3D62;
    padding: 4px;
}

QPushButton {
    background-color: #3498DB;
    color: #FFFFFF;
    border: none;
    border-radius: 8px;
    padding: 12px 24px;
    font-weight: 600;
    font-size: 14px;
    min-height: 20px;
}

QPushButton:hover {
    background-color: #2E86C1;
}

QPushButton:pressed {
    background-color: #21618C;
}

QPushButton#generateButton {
    background-color: #1B4F72;
    font-size: 16px;
    padding: 16px 32px;
    font-weight: bold;
    letter-spacing: 1px;
}

QPushButton#generateButton:hover {
    background-color: #154360;
}

QScrollArea {
    border: none;
    background-color: transparent;
}

QFrame#headerFrame {
    background-color: #FFFFFF;
    border-bottom: 3px solid #3498DB;
}

QFrame#sectionFrame {
    background-color: #FFFFFF;
    border: 2px solid #AED6F1;
    border-radius: 10px;
    margin: 8px 0px;
    padding: 12px;
}

QLabel#titleLabel {
    color: #0A3D62;
    font-size: 28px;
    font-weight: bold;
    letter-spacing: 1px;
}

QLabel#subtitleLabel {
    color: #5DADE2;
    font-size: 20px;
    font-weight: 500;
}

QLabel#sectionLabel {
    color: #1B4F72;
    font-size: 15px;
    font-weight: bold;
    padding: 8px 0px;
}


QLabel#fileLabel {
    color: #3498DB;
    background-color: #F4F9FD;
    border: 1px solid #AED6F1;
    border-radius: 6px;
    padding: 10px;
    font-size: 12px;
}
"""

# ---------------------------------------
#  MAIN GUI CLASS — MODERN VERSION
# ---------------------------------------

class GUI(QWidget):

    def __init__(self, root, user_settings):
        # Create QApplication BEFORE QWidget
        self._app = QApplication.instance()
        if self._app is None:
            self._app = QApplication(sys.argv)

        super().__init__()

        self.root = root
        self.user_settings = user_settings
        self.quit_status = False

        # Apply modern stylesheet
        self.setStyleSheet(MODERN_STYLESHEET)

        # Window settings
        self.setWindowTitle("DNABOT Configuration")
        self.setMinimumWidth(750)
        self.setMinimumHeight(700)





        # # Main layout
        # main_layout = QVBoxLayout(self)
        # main_layout.setContentsMargins(0, 0, 0, 0)
        # main_layout.setSpacing(0)

        # # Scroll area for content
        # scroll = QScrollArea()
        # scroll.setWidgetResizable(True)
        # scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        

        # # Header with logo
        # header = self.create_header()
        # main_layout.addWidget(header)


        # container = QWidget()
        # self.layout = QVBoxLayout(container)
        # self.layout.setContentsMargins(20, 20, 20, 20)
        # self.layout.setSpacing(15)
        
        # scroll.setWidget(container)
        # main_layout.addWidget(scroll)


        # Main layout
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)

        # Scroll area for content
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        container = QWidget()
        self.layout = QVBoxLayout(container)
        self.layout.setContentsMargins(20, 20, 20, 20)
        self.layout.setSpacing(15)

        # Header with logo (first thing in scrollable area)
        header = self.create_header()
        self.layout.addWidget(header)
        self.layout.addSpacing(10)

        scroll.setWidget(container)
        main_layout.addWidget(scroll)





        # Robot Type Section
        self.add_section_header("Robot Configuration")
        robot_frame = self.create_frame()
        robot_layout = QVBoxLayout(robot_frame)
        
        self.robot_type = self.add_dropdown_to_layout(
            robot_layout,
            "Robot Type",
            ["OT-2", "Flex"],
            callback=self.update_robot_type
        )
        
        self.layout.addWidget(robot_frame)

        # Hardware Section
        self.add_section_header("Hardware Setup")
        hardware_frame = self.create_frame()
        hardware_layout = QVBoxLayout(hardware_frame)
        
        self.single_pipette = self.add_dropdown_to_layout(
            hardware_layout, "Single Pipette", OT2_SINGLE_PIPETTES
        )
        self.multi_pipette = self.add_dropdown_to_layout(
            hardware_layout, "Multi Pipette", OT2_MULTI_PIPETTES
        )
        self.single_mount = self.add_dropdown_to_layout(
            hardware_layout, "Single Pipette Mount", ["left", "right"]
        )
        self.multi_mount = self.add_dropdown_to_layout(
            hardware_layout, "Multi Pipette Mount", ["left", "right"]
        )
        self.thermocycler = self.add_dropdown_to_layout(
            hardware_layout, "Thermocycler", ["thermocyclerModuleV1", "thermocyclerModuleV2"]
        )
        self.magdeck = self.add_dropdown_to_layout(
            hardware_layout, "Magnetic Module", 
            ["magnetic module gen1", "magnetic module gen2", "magneticBlockV1"]
        )
        
        self.layout.addWidget(hardware_frame)

        # Clip Parameters Section
        self.add_section_header("Clip Parameters")
        params_frame = self.create_frame()
        params_layout = QVBoxLayout(params_frame)

        params = self.user_settings["parameters"]
        self.premix_linkers = self.add_dropdown_to_layout(
            params_layout,
            "Premix Linkers",
            ["Yes", "No"]
        )
        self.premix_linkers.setCurrentText(str(params["premix_linkers"]["value"]))

        self.premix_parts = self.add_dropdown_to_layout(
            params_layout,
            "Premix Parts",
            ["Yes", "No"]
        )
        self.premix_parts.setCurrentText(str(params["premix_parts"]["value"]))

        self.linkers_volume = self.add_line_edit_to_layout(
            params_layout,
            "Linkers Volume (uL)",
            str(params["linkers_volume"]["value"])
        )

        self.parts_volume = self.add_line_edit_to_layout(
            params_layout,
            "Parts Volume (uL)",
            str(params["parts_volume"]["value"])
        )

        self.thermo_temp = self.add_line_edit_to_layout(
            params_layout,
            "Thermocycler Temp (C)",
            str(params["thermo_temp"]["value"])
        )

        self.layout.addWidget(params_frame)

        # Labware Section
        self.add_section_header("Labware Configuration")
        labware_frame = self.create_frame()
        labware_layout = QVBoxLayout(labware_frame)
        
        self.labware_widgets = {}
        for key in ALL_LABWARE_KEYS:
            items = OT2_LABWARE.get(key, FLEX_LABWARE.get(key, []))
            widget = self.add_dropdown_to_layout(labware_layout, key.replace("_", " ").title(), items)
            self.labware_widgets[key] = widget
        
        self.layout.addWidget(labware_frame)

        # File Selection Section
        self.add_section_header("File Selection")
        file_frame = self.create_frame()
        file_layout = QVBoxLayout(file_frame)
        
        self.construct_path = ""
        self.source_paths = []

        btn_construct = QPushButton("Select Construct CSV File")
        btn_construct.clicked.connect(self.select_construct)
        file_layout.addWidget(btn_construct)

        self.construct_label = QLabel("No construct file selected")
        self.construct_label.setObjectName("fileLabel")
        self.construct_label.setWordWrap(True)
        file_layout.addWidget(self.construct_label)

        btn_sources = QPushButton("Select BASIC Part CSV Files")
        btn_sources.clicked.connect(self.select_sources)
        file_layout.addWidget(btn_sources)

        self.sources_label = QLabel("No BASIC part files selected")
        self.sources_label.setObjectName("fileLabel")
        self.sources_label.setWordWrap(True)
        self.sources_label.setMinimumHeight(50)
        self.sources_label.setMaximumHeight(150)

        self.sources_scroll = QScrollArea()
        self.sources_scroll.setWidgetResizable(True)
        self.sources_scroll.setWidget(self.sources_label)
        self.sources_scroll.setMaximumHeight(150)
        file_layout.addWidget(self.sources_scroll)
        
        self.layout.addWidget(file_frame)

        # Generate button
        self.layout.addSpacing(20)
        generate_btn = QPushButton("GENERATE PROTOCOL")
        generate_btn.setObjectName("generateButton")
        generate_btn.clicked.connect(self.generate)
        self.layout.addWidget(generate_btn)
        generate_btn.setMaximumWidth(400)  # Add this line
        generate_btn.setMinimumWidth(300)  # Add this line

        # Center the button, Added this block
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(generate_btn)
        button_layout.addStretch()
        self.layout.addLayout(button_layout)

        self.layout.addStretch()

        self.update_robot_type()
        self.show()
        self._app.exec()

    def create_header(self):
        """Create the header with logo and title"""
        header_frame = QFrame()
        header_frame.setObjectName("headerFrame")
        header_frame.setFixedHeight(200)
        
        header_layout = QVBoxLayout(header_frame)  # Add, Changed to VBoxLayout
        header_layout.setContentsMargins(30, 20, 30, 20)
        header_layout.setSpacing(10)
        
        # Add stretch before to center content
        #header_layout.addStretch()
        
        # Logo
        logo_label = QLabel()
        logo_path = "images/dnabot_logo.png"
        
        pixmap = QPixmap(logo_path)
        scaled_pixmap = pixmap.scaled(200, 200, Qt.KeepAspectRatio, Qt.SmoothTransformation)
        
        logo_label.setPixmap(scaled_pixmap)
        
        logo_label.setAlignment(Qt.AlignCenter)
        header_layout.addWidget(logo_label)
        
        # Add spacing between logo and title
        #header_layout.addSpacing(20)
        
        # Title section
        # title_layout = QVBoxLayout()
        # title_layout.setSpacing(5)
        
        # title = QLabel("DNABot Configuration")
        # title.setObjectName("titleLabel")
        # title.setAlignment(Qt.AlignVCenter)
        # title_layout.addWidget(title)
        
        # subtitle = QLabel("Automated DNA Assembly Protocol Generator")
        # subtitle.setObjectName("subtitleLabel")
        # subtitle.setAlignment(Qt.AlignVCenter)
        # title_layout.addWidget(subtitle)
        
        # header_layout.addLayout(title_layout)
        
        # # Add stretch after to center content
        # header_layout.addStretch()

        # Title
        title = QLabel("Automated BASIC DNA Assembly Protocols")
        title.setObjectName("titleLabel")
        title.setAlignment(Qt.AlignCenter)
        header_layout.addWidget(title)
        
        # Subtitle
        subtitle = QLabel("Opentron Hackathon 2025")
        subtitle.setObjectName("subtitleLabel")
        subtitle.setAlignment(Qt.AlignCenter)
        header_layout.addWidget(subtitle)
        
        return header_frame

    def create_frame(self):
        """Create a styled frame for sections"""
        frame = QFrame()
        frame.setObjectName("sectionFrame")
        return frame

    def add_section_header(self, text):
        """Add a section header label"""
        label = QLabel(text)
        label.setObjectName("sectionLabel")
        self.layout.addWidget(label)

    def add_dropdown(self, label, options, callback=None):
        """Add a dropdown to the main layout"""
        form = QFormLayout()
        form.setSpacing(8)
        form.setContentsMargins(10, 5, 10, 5)
        
        label_widget = QLabel(label)
        label_widget.setStyleSheet("font-weight: 600;")
        
        widget = QComboBox()
        widget.addItems(options)
        if callback:
            widget.currentTextChanged.connect(callback)
        
        form.addRow(label_widget, widget)
        self.layout.addLayout(form)
        return widget

    def add_dropdown_to_layout(self, parent_layout, label, options, callback=None):
        """Add a dropdown to a specific layout"""
        form = QFormLayout()
        form.setSpacing(8)
        form.setContentsMargins(10, 5, 10, 5)
        
        label_widget = QLabel(label)
        label_widget.setStyleSheet("font-weight: 600;")
        
        widget = QComboBox()
        widget.addItems(options)
        if callback:
            widget.currentTextChanged.connect(callback)

        widget._label_widget = label_widget
        form.addRow(label_widget, widget)
        parent_layout.addLayout(form)
        return widget

    def add_line_edit_to_layout(self, parent_layout, label, value):
        """Add a line edit to a specific layout"""
        form = QFormLayout()
        form.setSpacing(8)
        form.setContentsMargins(10, 5, 10, 5)

        label_widget = QLabel(label)
        label_widget.setStyleSheet("font-weight: 600;")

        widget = QLineEdit()
        widget.setText(value)

        form.addRow(label_widget, widget)
        parent_layout.addLayout(form)
        return widget

    def update_robot_type(self):
        """Update hardware and labware defaults based on robot type"""
        robot = self.robot_type.currentText()
        if robot == "OT-2":
            lab = OT2_LABWARE
            defaults = OT2_DEFAULTS
            self.single_pipette.clear()
            self.single_pipette.addItems(OT2_SINGLE_PIPETTES)
            self.single_pipette.setCurrentIndex(0)
            self.multi_pipette.clear()
            self.multi_pipette.addItems(OT2_MULTI_PIPETTES)
            self.multi_pipette.setCurrentIndex(0)
        else:
            lab = FLEX_LABWARE
            defaults = FLEX_DEFAULTS
            self.single_pipette.clear()
            self.single_pipette.addItems(FLEX_SINGLE_PIPETTES)
            self.single_pipette.setCurrentIndex(0)
            self.multi_pipette.clear()
            self.multi_pipette.addItems(FLEX_MULTI_PIPETTES)
            self.multi_pipette.setCurrentIndex(0)

        self.single_pipette.setCurrentText(defaults["single_pipette"])
        self.multi_pipette.setCurrentText(defaults["multi_pipette"])
        self.single_mount.setCurrentText(defaults["single_mount"])
        self.multi_mount.setCurrentText(defaults["multi_mount"])
        self.thermocycler.setCurrentText(defaults["thermocycler"])
        self.magdeck.setCurrentText(defaults["magdeck"])

        visible_keys = set(lab.keys())
        for key, widget in self.labware_widgets.items():
            is_visible = key in visible_keys
            widget.setVisible(is_visible)
            widget._label_widget.setVisible(is_visible)
            if is_visible:
                widget.clear()
                widget.addItems(lab[key])
                widget.setCurrentIndex(0)

    def select_construct(self):
        """Open file dialog for construct CSV"""
        file, _ = QFileDialog.getOpenFileName(self, "Select Construct CSV")
        if file:
            self.construct_path = file
            filename = os.path.basename(file)
            self.construct_label.setText(f"✓ Selected: {filename}\n{file}")
        else:
            self.construct_label.setText("No construct file selected")

    def select_sources(self):
        """Open file dialog for source CSV files"""
        files, _ = QFileDialog.getOpenFileNames(self, "Select BASIC Part CSV Files")
        if files:
            self.source_paths = files
            file_list = "<b>✓ Selected files:</b><br><ul style='margin-top: 5px;'>"
            for file in files:
                filename = os.path.basename(file)
                file_list += f"<li>{filename}</li>"
            file_list += "</ul>"
            self.sources_label.setText(file_list)
        else:
            self.sources_label.setText("No BASIC part files selected")

    def generate(self):
        """Save settings and close"""
        us = self.user_settings

        robot = self.robot_type.currentText()
        us["robot_type"] = robot
        us["hardware"]["robot_type"]["id"] = robot

        us["hardware"]["single_pipette"]["id"] = self.single_pipette.currentText()
        us["hardware"]["multi_pipette"]["id"] = self.multi_pipette.currentText()
        us["hardware"]["thermocycler"]["id"] = self.thermocycler.currentText()
        us["hardware"]["mag_deck"]["id"] = self.magdeck.currentText()

        us["hardware"]["single_pipette_mount"]["id"] = self.single_mount.currentText()
        us["hardware"]["multi_pipette_mount"]["id"] = self.multi_mount.currentText()

        us["parameters"]["premix_linkers"]["value"] = self.premix_linkers.currentText()
        us["parameters"]["premix_parts"]["value"] = self.premix_parts.currentText()
        us["parameters"]["linkers_volume"]["value"] = float(self.linkers_volume.text())
        us["parameters"]["parts_volume"]["value"] = float(self.parts_volume.text())
        us["parameters"]["thermo_temp"]["value"] = float(self.thermo_temp.text())

        for key, widget in self.labware_widgets.items():
            us["labwares"][key]["id"] = widget.currentText()

        us["construct_path"] = self.construct_path
        us["sources_paths"] = self.source_paths

        self.quit_status = False
        self.close()
