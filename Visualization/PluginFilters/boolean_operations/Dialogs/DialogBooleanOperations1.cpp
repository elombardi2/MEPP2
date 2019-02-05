// Copyright (c) 2012-2019 University of Lyon and CNRS (France).
// All rights reserved.
//
// This file is part of MEPP2; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published 
// by the Free Software Foundation; either version 3 of the License, 
// or (at your option) any later version.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
#include "DialogBooleanOperations1.h"
#include "ui_DialogBooleanOperations1.h"
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogBooleanOperations1::DialogBooleanOperations1(QWidget *parent)
    : QDialog(parent), ui(new Ui::DialogBooleanOperations1)
{
  ui->setupUi(this);
}
////////////////////////////////////////////////////////////////////////////////
FEVV::DialogBooleanOperations1::~DialogBooleanOperations1() { delete ui; }
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogBooleanOperations1::setParameters(double x, double y, double z)
{
  ui->lineEdit_X->setText(QString::number(x));
  ui->lineEdit_Y->setText(QString::number(y));
  ui->lineEdit_Z->setText(QString::number(z));
}
////////////////////////////////////////////////////////////////////////////////
void
FEVV::DialogBooleanOperations1::getParameters(double &x, double &y, double &z)
{
  x = ui->lineEdit_X->text().toDouble();
  y = ui->lineEdit_Y->text().toDouble();
  z = ui->lineEdit_Z->text().toDouble();
}
////////////////////////////////////////////////////////////////////////////////
